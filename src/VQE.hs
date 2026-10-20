-- | Discrete noiseless VQE evaluated through VHamilton.
--
-- VHamilton itself needs no VQE-specific syntax: a variational process is an
-- ordinary 'Proc' containing parameterized bare actions.  VQE adds only
-- an external sequential ansatz description telling us which VHamilton
-- processes are executed one after another.
--
-- For a parameter precision r and integer tuple q, the valuation is
--
--   rho(alpha_j) = q_j / 2^r.
--
-- Each open VHamilton block is compiled under rho using the ordinary
-- Hamilton Trotter scheduling/halving mechanism.  Thus variation selects a member of a
-- process family; Trotter refinement controls only how that selected member is
-- compiled.
module VQE
  ( VBlock(..)
  , VAnsatz(..)
  , VQEResult(..)
  , valuationForParams
  , compileAnsatz
  , energyForCircuit
  , energyForParams
  , gridSearch
  ) where

import Numeric.Natural (Natural)
import Control.Monad (foldM, when)
import qualified Data.Aeson as A
import qualified Data.ByteString.Lazy.Char8 as BL
import qualified Data.Map.Strict as M
import qualified Data.Set as SSet
import System.Exit (ExitCode(..))
import System.IO (hFlush, stdout)
import System.Process (readProcessWithExitCode)
import Text.Printf (printf)

import Syntax (Proc, QCircuit, Valuation, freeParams)
import Compiler (compileTrotterWith)
import Qasm (emit)

--------------------------------------------------------------------------------
-- VHamilton ansatz blocks

-- | One sequential ansatz block.  'vbProcess' may contain any number of free
-- VHamilton parameters; the block itself is still just a core Proc.
data VBlock = VBlock
  { vbProcess    :: Proc
  , vbRefinement :: Natural
  } deriving (Eq, Show)

-- | External VQE metadata for a family of VHamilton processes.
--
-- Parameter precision is deliberately not syntax: the same VHamilton process
-- family can be searched at r=3, r=7, ... without changing the term.
data VAnsatz = VAnsatz
  { vaParamNames :: [String]
  , vaPrecision  :: Natural
  , vaBlocks     :: [VBlock]
  } deriving (Eq, Show)

-- | Convert the discrete integer search coordinates to the real/dyadic
-- valuation used by VHamilton: rho(alpha_j) = q_j / 2^r.
valuationForParams :: VAnsatz -> [Integer] -> Valuation
valuationForParams ansatz qs
  | length names /= length qs =
      error "valuationForParams: number of parameter names and values must match"
  | length names /= SSet.size (SSet.fromList names) =
      error "valuationForParams: duplicate parameter names"
  | SSet.fromList names /= referenced =
      error ("valuationForParams: declared parameters " ++ show names ++
             " differ from parameters in the VHamilton blocks " ++ show (SSet.toAscList referenced))
  | otherwise = M.fromList (zip names values)
  where
    names = vaParamNames ansatz
    referenced = SSet.fromList (concatMap (freeParams . vbProcess) (vaBlocks ansatz))
    denom = fromIntegral
      (2 ^ (fromIntegral (vaPrecision ansatz) :: Int) :: Integer)
    values = map ((/ denom) . fromIntegral) qs

-- | Compile a sequential VQE ansatz. Each open block is interpreted under
-- the supplied valuation and compiled at its own Trotter refinement level;
-- no arbitrary real literal is inserted into the object language.
compileAnsatz :: VAnsatz -> [Integer] -> QCircuit
compileAnsatz ansatz qs =
  mconcat
    [ compileTrotterWith rho (vbRefinement block) (vbProcess block)
    | block <- vaBlocks ansatz
    ]
  where
    rho = valuationForParams ansatz qs

--------------------------------------------------------------------------------
-- Energy oracle (delegated to Python/Qiskit, like the existing Denot module)

data PyEnergy = PyEnergy
  { peEnergy :: Double
  , peGround :: Double
  }

instance A.FromJSON PyEnergy where
  parseJSON = A.withObject "PyEnergy" $ \o ->
    PyEnergy <$> o A..: "energy" <*> o A..: "ground_energy"

-- | Exact statevector energy of a compiled variational circuit against a
-- closed Hamilton target Proc.  Also returns the exact ground-state energy of
-- the target, useful for validating the small examples.
energyForCircuit :: Proc -> QCircuit -> IO (Double, Double)
energyForCircuit target qc = do
  qasm <- emit qc
  let payload = BL.unpack $ A.encode $ A.object
        [ "hamiltonian" A..= target
        , "qasm"        A..= qasm
        ]
  (ec, out, err) <-
    readProcessWithExitCode "python3" ["-m", "scripts.vqe.vqe_eval"] payload
  case ec of
    ExitFailure k ->
      error ("scripts.vqe.vqe_eval failed: " ++ show k ++ "\nSTDERR:\n" ++ err)
    ExitSuccess ->
      case A.eitherDecode (BL.pack out) of
        Left e ->
          error ("scripts.vqe.vqe_eval bad JSON: " ++ e
                 ++ "\nSTDERR:\n" ++ err
                 ++ "\nSTDOUT:\n" ++ out)
        Right r -> pure (peEnergy r, peGround r)

--------------------------------------------------------------------------------
-- Classical discrete optimization

data VQEResult = VQEResult
  { vrParams       :: [Integer]
  , vrEnergy       :: Double
  , vrGroundEnergy :: Double
  } deriving (Eq, Show)

energyForParams :: Proc -> VAnsatz -> [Integer] -> IO VQEResult
energyForParams target ansatz qs = do
  let qc = compileAnsatz ansatz qs
  (e, e0) <- energyForCircuit target qc
  pure VQEResult
    { vrParams = qs
    , vrEnergy = e
    , vrGroundEnergy = e0
    }

-- | Exhaustive optimizer over a finite integer grid.  The grid is classical;
-- every energy query is produced by instantiating and compiling VHamilton.
gridSearch
  :: Proc
  -> VAnsatz
  -> [[Integer]]
  -> IO (VQEResult, [VQEResult])
gridSearch target ansatz grids = do
  let tuples = sequence grids
      total = length tuples
      interval = max 1 (total `div` 100)

      step
        :: (Int, [VQEResult], Maybe VQEResult)
        -> [Integer]
        -> IO (Int, [VQEResult], Maybe VQEResult)
      step (i0, acc, maybeBest) qs = do
        result <- energyForParams target ansatz qs
        let i = i0 + 1
            best =
              case maybeBest of
                Nothing -> result
                Just old
                  | vrEnergy result < vrEnergy old -> result
                  | otherwise                      -> old
            pct = 100.0 * fromIntegral i / fromIntegral total :: Double

        when (i == 1 || i `mod` interval == 0 || i == total) $ do
          printf
            "\rVQE progress: %d/%d (%5.1f%%)  current q=%s  E=%.12f  best q=%s  best E=%.12f"
            i
            total
            pct
            (show (vrParams result))
            (vrEnergy result)
            (show (vrParams best))
            (vrEnergy best)
          hFlush stdout

        pure (i, result : acc, Just best)

  if null tuples
    then error "gridSearch: empty parameter grid"
    else do
      (_, reversedResults, maybeBest) <- foldM step (0, [], Nothing) tuples
      putStrLn ""
      case maybeBest of
        Nothing -> error "gridSearch: no VQE results"
        Just best -> pure (best, reverse reversedResults)
