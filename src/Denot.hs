-- | Verification oracle: compare compiled circuits to ideal unitaries.
--
-- Delegates to scripts.semantics.denot, which implements the Hamiltonian and unitary
-- semantics:
--   ⟦p⟧_H (hamiltonian) and ⟦p⟧_U = exp(-2πi ⟦p⟧_H) (ideal unitary).
--
-- compareToIdeal   -- ‖(⟦sch(p)/2^n⟧_U)^{2^n} - ⟦p⟧_U‖_F  (Trotter error)
-- compareIdeals    -- ‖⟦p1⟧_U - ⟦p2⟧_U‖_F  (dyadic approximation error)
module Denot where

import GHC.Generics (Generic)
import qualified Data.Aeson as A
import qualified Data.ByteString.Lazy.Char8 as BL
import System.Exit (ExitCode(..))
import System.Process (readProcessWithExitCode)

import Syntax (Proc, QCircuit, Valuation)
import qualified Data.Map.Strict as M
import Qasm (emit)

newtype IdealFrobResult = IdealFrobResult { idealFrob :: Double }

instance A.FromJSON IdealFrobResult where
  parseJSON = A.withObject "IdealFrobResult" $ \o ->
    IdealFrobResult <$> o A..: "frob"

data DenotResult = DenotResult
  { frobenius :: Double
  , idealU    :: [[(Double, Double)]]
  } deriving (Eq, Show, Generic)

instance A.FromJSON DenotResult where
  parseJSON = A.withObject "DenotResult" $ \o ->
    DenotResult <$> o A..: "frob_aligned" <*> o A..: "ideal"

-- | Run scripts.semantics.denot on a (process, circuit) pair.
-- Returns the phase-aligned Frobenius distance and the ideal matrix.
compareToIdeal
  :: Proc
  -> QCircuit
  -> IO (Double, [[(Double, Double)]])
compareToIdeal = compareToIdealWith M.empty

-- | Compare an open VHamilton process to its compiled circuit at valuation rho.
compareToIdealWith :: Valuation -> Proc -> QCircuit -> IO (Double, [[(Double, Double)]])
compareToIdealWith rho proc qc = do
  qasm <- emit qc
  let payload = BL.unpack $ A.encode $ A.object
        [ "proc" A..= proc
        , "valuation" A..= rho
        , "qasm" A..= qasm
        ]
  (ec, out, err) <-
    readProcessWithExitCode "python3" ["-m", "scripts.semantics.denot"] payload
  case ec of
    ExitFailure k ->
      error ("scripts.semantics.denot failed: " ++ show k ++ "\nSTDERR:\n" ++ err)
    ExitSuccess ->
      case A.eitherDecode (BL.pack out) of
        Left e ->
          error ("scripts.semantics.denot bad JSON: " ++ e
                 ++ "\nSTDERR:\n" ++ err
                 ++ "\nSTDOUT:\n" ++ out)
        Right r -> pure (frobenius r, idealU r)

-- | Compare two processes' ideal unitaries without approximation.
-- Returns ‖exp(-2πi ⟦p1⟧) - exp(-2πi ⟦p2⟧)‖_F (no phase alignment).
compareIdeals :: Proc -> Proc -> IO Double
compareIdeals p1 p2 = do
  let payload = BL.unpack $ A.encode $ A.object
        [ "proc1" A..= p1
        , "proc2" A..= p2
        ]
  (ec, out, err) <-
    readProcessWithExitCode "python3" ["-m", "scripts.semantics.compare_ideals"] payload
  case ec of
    ExitFailure k ->
      error ("scripts.semantics.compare_ideals failed: " ++ show k ++ "\nSTDERR:\n" ++ err)
    ExitSuccess ->
      case A.eitherDecode (BL.pack out) of
        Left e  ->
          error ("scripts.semantics.compare_ideals bad JSON: " ++ e
                 ++ "\nSTDERR:\n" ++ err ++ "\nSTDOUT:\n" ++ out)
        Right r -> pure (idealFrob r)

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
