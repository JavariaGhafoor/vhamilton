-- | Regression tests for the Hamilton compiler and the first discrete VQE demo.
--
-- These tests deliberately exercise non-self-adjoint/order-sensitive Clifford
-- structure (especially Y via S · H) and two-qubit tensor placement.  They are
-- intended to catch the compiler ordering bug that previously turned +Y into
-- -Y, as well as future qubit-order regressions.
module Main where

import Control.Monad (forM_)
import qualified Data.Map.Strict as M
import qualified Data.Aeson as A
import qualified Data.ByteString.Lazy.Char8 as BL
import System.Exit (ExitCode(..), exitFailure)
import System.Process (readProcessWithExitCode)
import Text.Printf (printf)

import Compiler (compileTrotter, compileTrotterWith, diagonalAnglesWith)
import Denot (compareToIdeal, compareToIdealWith)
import ExamplesChem (ChemTerm(..), chemTerm)
import Syntax
  ( Clifford(..), (·)
  , Proc(..), (∥)
  , BareAction(..), deltaPow, (⊕), Action(..), freeParams
  )
import VQE (VBlock(..), VAnsatz(..), VQEResult(..), gridSearch)
import Semantics (seg0)

--------------------------------------------------------------------------------
-- Small constructors

act :: BareAction -> Proc
act = Atom . Bare

conj :: Clifford -> BareAction -> Proc
conj u d = Atom (Conj u (Bare d))

-- Full Pauli generators (eigenvalues +/-1), used by the VQE target/generator.
pauliZ :: Proc
pauliZ = act (deltaPow 0 1 ⊕ deltaPow 0 (-1))

pauliX :: Proc
pauliX = conj H (deltaPow 0 1 ⊕ deltaPow 0 (-1))

pauliY :: Proc
pauliY = conj (S · H) (deltaPow 0 1 ⊕ deltaPow 0 (-1))

-- Weighted Pauli generators at 1/8.  Using a non-integral dyadic coefficient is
-- important: exp(-2*pi*i*P) would be identity for eigenvalues +/-1 and would
-- therefore be too weak a compiler regression test.
weightedZ :: Proc
weightedZ = act (deltaPow 3 1 ⊕ deltaPow 3 (-1))

weightedX :: Proc
weightedX = conj H (deltaPow 3 1 ⊕ deltaPow 3 (-1))

weightedY :: Proc
weightedY = conj (S · H) (deltaPow 3 1 ⊕ deltaPow 3 (-1))

weighted2 :: String -> Proc
weighted2 lbl = chemTerm (ChemTerm lbl 1 3)

--------------------------------------------------------------------------------
-- Assertions

tol :: Double
tol = 1e-9

failTest :: String -> String -> IO a
failTest name msg = do
  putStrLn $ "FAIL  " ++ name ++ ": " ++ msg
  exitFailure

assertClose :: String -> Double -> Double -> IO ()
assertClose name actual expected =
  if abs (actual - expected) <= tol
    then putStrLn $ "PASS  " ++ name
    else failTest name
      ("expected " ++ printf "%.12g" expected
       ++ ", got " ++ printf "%.12g" actual)

assertEqual :: (Eq a, Show a) => String -> a -> a -> IO ()
assertEqual name actual expected =
  if actual == expected
    then putStrLn $ "PASS  " ++ name
    else failTest name
      ("expected " ++ show expected ++ ", got " ++ show actual)

--------------------------------------------------------------------------------
-- Compiler-vs-semantics tests

-- Check U_compiled == exp(-2*pi*i [[p]]_H), up to global phase.
checkCompiled :: String -> Proc -> IO ()
checkCompiled name p = do
  let qc = compileTrotter 0 p
  (frob, _) <- compareToIdeal p qc
  if frob <= tol
    then putStrLn $ "PASS  compiled " ++ name
                    ++ "  (phase-aligned Frobenius=" ++ printf "%.3e" frob ++ ")"
    else failTest ("compiled " ++ name)
      ("phase-aligned Frobenius distance " ++ printf "%.12g" frob)

-- Check that the declared chemistry label itself denotes the standard matrix
-- coefficient * (P_first tensor P_second), independently of compiled QASM.
-- This makes XI/IX and XY/YX wire placement explicit rather than merely
-- checking internal consistency between two conventions.
data LabelResult = LabelResult { labelFrob :: Double }

instance A.FromJSON LabelResult where
  parseJSON = A.withObject "LabelResult" $ \o ->
    LabelResult <$> o A..: "frob"

checkLabel :: String -> Proc -> IO ()
checkLabel lbl p = do
  let payload = BL.unpack $ A.encode $ A.object
        [ "proc" A..= p
        , "pauli" A..= lbl
        , "coeff" A..= (0.125 :: Double)
        ]
  (ec, out, err) <-
    readProcessWithExitCode "python3" ["-m", "scripts.semantics.pauli_semantics_check"] payload
  case ec of
    ExitFailure k -> failTest ("label " ++ lbl)
      ("scripts.semantics.pauli_semantics_check failed: " ++ show k ++ "\n" ++ err)
    ExitSuccess ->
      case A.eitherDecode (BL.pack out) of
        Left e -> failTest ("label " ++ lbl)
          ("bad JSON: " ++ e ++ "\nSTDOUT:\n" ++ out ++ "\nSTDERR:\n" ++ err)
        Right r ->
          if labelFrob r <= tol
            then putStrLn $ "PASS  label " ++ lbl
                            ++ "  (Hamiltonian Frobenius="
                            ++ printf "%.3e" (labelFrob r) ++ ")"
            else failTest ("label " ++ lbl)
              ("Hamiltonian Frobenius distance " ++ printf "%.12g" (labelFrob r))

--------------------------------------------------------------------------------
-- Minimal VHamilton syntax, halving, and explicit valuation checks

testMinimalSyntax :: IO ()
testMinimalSyntax = do
  let rho = M.fromList [("alpha", 3/16)]
      core = DeltaParam 1 "alpha" ⊕ NegDeltaParam 1 "alpha"
      p = conj (S · H) core
  assertEqual "parameter names" (freeParams p) ["alpha"]
  assertEqual "halving positive action" (seg0 (DeltaParam 1 "alpha") 2)
    (DeltaParam 3 "alpha")
  assertEqual "halving negative action" (seg0 (NegDeltaParam 1 "alpha") 2)
    (NegDeltaParam 3 "alpha")
  assertClose "positive parameter phase" (head (diagonalAnglesWith rho (DeltaParam 1 "alpha")))
    (-2 * pi * (3/16) / 2)
  assertClose "negative parameter phase" (head (diagonalAnglesWith rho (NegDeltaParam 1 "alpha")))
    (2 * pi * (3/16) / 2)
  let qc = compileTrotterWith rho 0 p
  (frob, _) <- compareToIdealWith rho p qc
  assertClose "open VHamilton compiled vs ideal" frob 0

--------------------------------------------------------------------------------
-- VQE regression

testVQE :: IO ()
testVQE = do
  let target = pauliX ∥ pauliZ
      pauliYVar =
        conj (S · H)
          (DeltaParam 0 "alpha" ⊕ NegDeltaParam 0 "alpha")
      ansatz = VAnsatz
        { vaParamNames = ["alpha"]
        , vaPrecision = 4
        , vaBlocks = [VBlock pauliYVar 0]
        }
  (best, _) <- gridSearch target ansatz [[0..7]]
  assertEqual "VQE best q" (vrParams best) [5]
  assertClose "VQE best energy" (vrEnergy best) (-(sqrt 2))
  assertClose "VQE exact ground" (vrGroundEnergy best) (-(sqrt 2))


--------------------------------------------------------------------------------

main :: IO ()
main = do
  putStrLn "Hamilton compiler/VQE regression tests"
  putStrLn ""

  -- One-qubit basis/conjugation tests.
  forM_
    [ ("X/8", weightedX)
    , ("Y/8", weightedY)
    , ("Z/8", weightedZ)
    ] (uncurry checkCompiled)

  -- Two-qubit compiled-unitary tests.
  let labels = ["XI", "IX", "XY", "YX", "ZZ"]
  forM_ labels $ \lbl -> checkCompiled (lbl ++ "/8") (weighted2 lbl)

  -- Independent semantic label/wire-order checks.
  forM_ labels $ \lbl -> checkLabel lbl (weighted2 lbl)

  -- Check the new atomic syntax and its interpretation.
  testMinimalSyntax

  -- Freeze the analytically solved toy VQE as a regression test.
  testVQE


  putStrLn ""
  putStrLn "All regression tests passed."
