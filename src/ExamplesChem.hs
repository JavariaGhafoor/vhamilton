-- | DMANH+ hydrogen-bond dynamics.
--
-- Implements the following construction:
--   diagCore   -- diagonal core of a weighted Pauli string (BareAction layer)
--   basisCliff -- Clifford rotating all-Z diagonal into Pauli basis (H for X, H;S for Y)
--   chemTerm   -- one Pauli term as a Hamilton Atom
--   chemProc   -- full DMANH+ Hamiltonian as a parallel composition
--
-- Input comes from scripts.chemistry.chem_decompose, which Pauli-decomposes
-- the published matrix and rounds each coefficient to a dyadic level.
module ExamplesChem where

-- | Chemistry application: the DMANH+ hydrogen-bond Hamiltonian.
--
-- The Hamiltonian lives in exactly one place: scripts.chemistry.chem_decompose, which holds
-- the published matrix (Richerme et al., J. Phys. Chem. Lett. 14, 7256, 2023), Pauli-decomposes
-- it, and rounds each coefficient to a chosen dyadic level. This module fetches
-- that decomposition at runtime (mirroring how Denot/Qasm shell out to their
-- Python helpers) and assembles the corresponding S_k process.
--
-- Each weighted Pauli string (num/2^k)·P is realized as a SINGLE diagonal
-- action using Hamilton's δ_k^j generator.  Open VHamilton terms instead
-- contain δ_k^alpha or -δ_k^alpha, with the sign an atomic syntactic choice. X/Y factors are handled by conjugating the all-Z diagonal into the
-- proper basis (H for X, S·H for Y). This avoids any blow-up: the schedule has
-- one action per Pauli term, independent of the size of the integer numerators.

import GHC.Generics (Generic)
import qualified Data.Aeson as A
import qualified Data.ByteString.Lazy.Char8 as BL
import System.Exit (ExitCode(..))
import System.Process (readProcessWithExitCode)

import Numeric.Natural (Natural)

import Syntax
  ( Clifford(..), (⊗), (·)
  , Proc(..), (∥)
  , BareAction(..), Action(..)
  )
import Examples (runAndCompare)

--------------------------------------------------------------------------------
-- Fetching the decomposition from scripts.chemistry.chem_decompose
--
-- Protocol (scripts.chemistry.chem_decompose --block <b> --level <k>, JSON on stdout):
--   { "scale": float, "level": int,
--     "identity_coeff_mh": float, "scaled_identity": float,
--     "identity_num": int, "identity_k": int,
--     "physical_ground_mh": float, "scaled_physical_ground": float,
--     "terms": [ {"pauli": "ZZ", "num": int, "k": int}, ... ] }

data ChemTerm = ChemTerm
  { ctPauli :: String
  , ctNum   :: Integer
  , ctK     :: Natural
  } deriving (Eq, Show, Generic)

instance A.FromJSON ChemTerm where
  parseJSON = A.withObject "ChemTerm" $ \o ->
    ChemTerm <$> o A..: "pauli"
             <*> o A..: "num"
             <*> (fromInteger <$> o A..: "k")

data ChemDecomp = ChemDecomp
  { cdScale                :: Double
  , cdLevel                :: Natural
  , cdIdentityCoeffMh      :: Double
  , cdScaledIdentity       :: Double
  , cdIdentityNum          :: Integer
  , cdIdentityK            :: Natural
  , cdPhysicalGroundMh     :: Double
  , cdScaledPhysicalGround :: Double
  , cdTerms                :: [ChemTerm]
  } deriving (Eq, Show, Generic)

instance A.FromJSON ChemDecomp where
  parseJSON = A.withObject "ChemDecomp" $ \o ->
    ChemDecomp <$> o A..: "scale"
               <*> (fromInteger <$> o A..: "level")
               <*> o A..: "identity_coeff_mh"
               <*> o A..: "scaled_identity"
               <*> o A..: "identity_num"
               <*> (fromInteger <$> o A..: "identity_k")
               <*> o A..: "physical_ground_mh"
               <*> o A..: "scaled_physical_ground"
               <*> o A..: "terms"

fetchDecomp :: String -> Natural -> IO ChemDecomp
fetchDecomp block level = do
  (ec, out, err) <- readProcessWithExitCode "python3"
    [ "-m", "scripts.chemistry.chem_decompose"
    , "--block", block
    , "--level", show (toInteger level)
    ] ""
  case ec of
    ExitFailure k ->
      error ("scripts.chemistry.chem_decompose failed: " ++ show k ++ "\n" ++ err)
    ExitSuccess ->
      case A.eitherDecode (BL.pack out) of
        Left e  -> error ("scripts.chemistry.chem_decompose bad JSON: " ++ e ++ "\n" ++ out)
        Right d -> pure d

--------------------------------------------------------------------------------
-- Building a weighted Pauli string as a single S_k action
--
-- For a 2-qubit label `lbl` (e.g. "ZI", "XX", "YY"), the diagonal of the
-- Z-projected Pauli string (X,Y -> Z) over basis states |00>,|01>,|10>,|11> has
-- sign exponent  par(i) = (number of Z-marked qubits in state |1>) mod 2.
-- Entry i carries Hamilton's integer dyadic action δ_k^j with
-- j = (-1)^par(i) · num. Symbolic entries use the explicit signed VHamilton
-- actions instead.

-- | Signs of the four diagonal entries of a two-qubit Pauli string after
-- projecting X/Y axes to Z.
diagSigns :: String -> [Double]
diagSigns lbl =
  [ if par i j == 0 then 1 else -1 | (i, j) <- basisIndices ]
  where
    basisIndices = [(0,0), (0,1), (1,0), (1,1)] :: [(Int,Int)]
    zmask = map (`elem` ("XYZ" :: String)) lbl
    par i j = sum [ b | (b, m) <- zip [i, j] zmask, m ] `mod` 2

-- | Build the four diagonal entries from a positive/negative choice.
-- The same sign pattern works for closed integer Hamilton generators and
-- symbolic VHamilton generators, without any scalar expressions in syntax.
diagCoreSigned :: String -> (Bool -> BareAction) -> BareAction
diagCoreSigned lbl makeEntry =
  case map (makeEntry . (> 0)) (diagSigns lbl) of
    [a, b, c, d] -> (a :⊕: b) :⊕: (c :⊕: d)
    _            -> error "diagCoreSigned: expected a 2-qubit (4-entry) label"

-- | Closed Hamilton dyadic core: integer numerator / 2^k.
diagCore :: Natural -> String -> Integer -> BareAction
diagCore k lbl num = diagCoreSigned lbl (\positive ->
  DeltaPow k (if positive then num else -num))

-- | VHamilton core: either +alpha or -alpha (with an optional overall sign).
diagCoreParam :: Natural -> String -> Bool -> String -> BareAction
diagCoreParam k lbl positiveAlpha alpha =
  diagCoreSigned lbl (\positive ->
    if positive == positiveAlpha then DeltaParam k alpha
                                 else NegDeltaParam k alpha)

-- | The Clifford rotating the all-Z diagonal into the basis of `lbl`.
--   X factor -> H, Y factor -> S·H, I/Z factor -> I.
--
--   NOTE on qubit order: signedExps places the diagonal sign pattern reading the
--   FIRST label character as the high (slow) bit and the SECOND as the low (fast)
--   bit. The pipeline's u ⊗ v applies u to the fast bit, so to land each Clifford
--   on the same qubit its phase pattern occupies, the tensor factors are built in
--   REVERSED label order: factor(q) ⊗ factor(p) for label "pq". (Verified against
--   xLow2 and the IX term: the forward order p ⊗ q crosses the wires and turns an
--   X term into a Z term.)
basisCliff :: String -> Maybe Clifford
basisCliff [p, q] =
  case (factor p, factor q) of
    (Nothing, Nothing) -> Nothing
    (mp, mq)           -> Just (orI mq ⊗ orI mp)   -- reversed: q on fast bit, p on slow
  where
    factor 'X' = Just H
    factor 'Y' = Just (S · H)
    factor _   = Nothing
    orI = maybe I id
basisCliff _ = Nothing

pauliCoreProc :: String -> BareAction -> Proc
pauliCoreProc lbl core =
  case basisCliff lbl of
    Nothing -> Atom (Bare core)
    Just u  -> Atom (Conj u (Bare core))

chemTerm :: ChemTerm -> Proc
chemTerm (ChemTerm lbl num k) = pauliCoreProc lbl (diagCore k lbl num)

chemProc :: ChemDecomp -> Proc
chemProc = foldr1 (∥) . map chemTerm . cdTerms

-- | Exact Pauli-string generator P with coefficient +1.
pauliProc :: String -> Proc
pauliProc lbl = pauliCoreProc lbl (diagCore 0 lbl 1)

-- | Symbolic +alpha/2^k times Pauli P, expressed only with VHamilton
-- parameterized bare actions.
pauliParamProc :: String -> Natural -> String -> Proc
pauliParamProc lbl k alpha = pauliSignedParamProc lbl k True alpha

-- | A negative coefficient is represented by choosing the opposite atomic
-- sign, not by constructing the forbidden scalar expression -alpha.
pauliSignedParamProc :: String -> Natural -> Bool -> String -> Proc
pauliSignedParamProc lbl k positive alpha =
  pauliCoreProc lbl (diagCoreParam k lbl positive alpha)

-- | The dyadic identity offset retained for VQE absolute-energy reporting.
-- This term is intentionally not compiled into chemProc because it would only
-- produce a global phase in circuit evolution.
identityDyadic :: ChemDecomp -> Double
identityDyadic d =
  fromIntegral (cdIdentityNum d) /
  fromIntegral (2 ^ (fromIntegral (cdIdentityK d) :: Int) :: Integer)

-- | Full dyadic chemistry Hamiltonian for energy-estimation/VQE semantics,
-- including the identity coefficient.  Unlike time evolution, VQE absolute
-- energies are shifted by c_I, so the identity term belongs in the target
-- Hamiltonian even though it need not be compiled into dynamics.
chemProcWithIdentity :: ChemDecomp -> Proc
chemProcWithIdentity d =
  chemTerm (ChemTerm "II" (cdIdentityNum d) (cdIdentityK d)) ∥ chemProc d

--------------------------------------------------------------------------------
-- Runners

-- | Fetch the upper block at coefficient level k, build the process, compile its
--   canonical Trotter circuit at refinement n, and report the Frobenius distance
--   to the exact dynamics.
testDMANH :: Natural -> Natural -> IO ()
testDMANH level n = do
  d <- fetchDecomp "upper" level
  runAndCompare False "chem.qasm" n (chemProc d)

-- | Sweep refinement levels at a fixed coefficient precision.
sweepDMANH :: Natural -> [Natural] -> IO ()
sweepDMANH level = mapM_ (testDMANH level)

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
