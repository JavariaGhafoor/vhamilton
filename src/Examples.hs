-- | Building blocks and small examples (§Examples, §Implementation).
--
-- Building blocks (§Examples):
--   zK / xK / yK   -- single-qubit Pauli generators at precision k; the basic
--                     Hamilton processes introduced in §Examples
-- Multi-term processes (§Examples, used in ordering/convergence tables):
--   zxK            -- Z ∥ X, asymmetric 2-term example
--   zxyK           -- Z ∥ X ∥ Y, 3-term example; 6 schedule orderings
--   zzParity       -- ZZ, commuting 2-qubit term (Trotter error = 0)
module Examples where

import Numeric.Natural (Natural)
import Data.List (sortOn)
import System.Exit (ExitCode(..))
import qualified Data.Text.IO as TIO
import System.Process (readProcessWithExitCode)
import Text.Printf (printf)

import Syntax
  ( Clifford(..), (⊗), (·)
  , Proc(..), (∥)
  , BareAction(..), delta, negDelta, (⊕), Action(..)
  , Schedule, QCircuit
  , prettyProc, prettySchedule
  )
import Compiler (compileTrotter, compileScheduleN)
import Semantics (canonicalSchedule, schedulePerms)
import Qasm (emitWithMeasurements)
import Denot (compareToIdeal)

--------------------------------------------------------------------------------
-- Formatting

fmtC :: (Double, Double) -> String
fmtC (re0, im0) =
  let eps     = 5e-3
      clamp x = if abs x < eps then 0 else x
      re      = clamp re0
      im      = clamp im0
      s       = case (re == 0, im == 0) of
                  (True,  True)  -> "0"
                  (False, True)  -> printf "%.2f" re
                  (True,  False) -> printf "%.2f" im ++ "i"
                  (False, False) ->
                    printf "%.2f" re ++
                    (if im < 0 then "" else "+") ++
                    printf "%.2f" im ++ "i"
  in printf "%12s" s

printMatrix :: [[(Double, Double)]] -> IO ()
printMatrix = mapM_ $ \row ->
  putStrLn $ "[" ++ concatMap (\z -> " " ++ fmtC z) row ++ " ]"

--------------------------------------------------------------------------------
-- Running circuits

-- | Emit, write, and optionally display a circuit via scripts.compiler.runq.
runCircuit :: Bool -> FilePath -> QCircuit -> IO ()
runCircuit showOutput fp qc = do
  qasmText <- emitWithMeasurements qc
  writeFile fp qasmText
  (ec, out, err) <-
    readProcessWithExitCode "python3" ["-m", "scripts.compiler.runq", fp] ""
  if showOutput then putStr out else pure ()
  if showOutput && not (null err) then putStr err else pure ()
  case ec of
    ExitSuccess   -> pure ()
    ExitFailure k -> error ("scripts.compiler.runq failed with exit code " ++ show k)

-- | Compile the Trotter circuit at level n, print ideal unitary and
-- Frobenius distance, and optionally display the circuit.
runAndCompare :: Bool -> FilePath -> Natural -> Proc -> IO ()
runAndCompare showCircuit fp n proc = do
  putStrLn $ "=== Proc ===\n" ++ prettyProc proc ++ "\n"
  putStrLn $ "=== Canonical schedule (level " ++ show n
          ++ ", repeated " ++ show (2 ^ n :: Int) ++ " times) ==="
  TIO.putStr (prettySchedule (canonicalSchedule proc n))
  let qc = compileTrotter n proc
  (frob, ui) <- compareToIdeal proc qc
  putStrLn "\n=== U_ideal ==="
  printMatrix ui
  putStrLn ""
  runCircuit showCircuit fp qc
  putStrLn $ "Frobenius: " ++ printf "%.4f" frob ++ "\n"

-- | Score one permutation against the ideal, optionally printing the circuit.
evalSchedule
  :: FilePath -> Natural -> Proc -> (Int, Schedule) -> IO (Double, Int)
evalSchedule fp n proc (i, sched) = do
  putStrLn $ "--- Permutation " ++ show i ++ " ---"
  TIO.putStr (prettySchedule sched)
  let qc = compileScheduleN n sched
  (frob, _) <- compareToIdeal proc qc
  runCircuit True fp qc
  putStrLn $ "Frobenius: " ++ printf "%.4f" frob ++ "\n"
  pure (frob, i)

-- | Run every permutation of the canonical schedule and rank by Frobenius.
comparePerms :: FilePath -> Natural -> Proc -> IO ()
comparePerms fp n proc = do
  putStrLn $ "=== Proc ===\n" ++ prettyProc proc ++ "\n"
  let perms = schedulePerms proc n
  putStrLn $ "Level: " ++ show n
          ++ "  |  Permutations: " ++ show (length perms) ++ "\n"
  results <- mapM (evalSchedule fp n proc) (zip [1..] perms)
  putStrLn "=== Ranking (best first) ==="
  mapM_ (\(frob, i) ->
    putStrLn $ "  #" ++ show (i :: Int)
            ++ "  Frobenius: " ++ printf "%.4f" frob)
    (sortOn fst results)
  putStrLn ""

--------------------------------------------------------------------------------
-- Smart constructors

act :: BareAction -> Proc
act = Atom . Bare

conj :: Clifford -> BareAction -> Proc
conj u d = Atom (Conj u (Bare d))

--------------------------------------------------------------------------------
-- Building blocks

-- §§ [building-blocks]
zK :: Natural -> Proc
zK k = act (Id ⊕ delta k)

xK :: Natural -> Proc
xK k = conj H (Id ⊕ delta k)

yK :: Natural -> Proc
yK k = conj (S · H) (Id ⊕ delta k)
-- §§ [/building-blocks]

xDagK :: Natural -> Proc
xDagK k = conj H (Id ⊕ negDelta k)

-- | 2-qubit diagonal with explicit phases on |00⟩,|01⟩,|10⟩,|11⟩.
--   k=0 entries use Id (angle 0, unaffected by segmentation);
--   k>0 entries use delta k (angle -2π/2^k).
diag2 :: Natural -> Natural -> Natural -> Natural -> Proc
diag2 a b c d = act ((f a ⊕ f b) ⊕ (f c ⊕ f d))
  where
    f 0 = Id
    f k = delta k

-- | Phase on |11⟩ only (CZ-like diagonal).
czLike :: Natural -> Proc
czLike k = diag2 0 0 0 k

-- §§ [zzParity]
zzParity :: Natural -> Proc
zzParity k = diag2 0 k k 0
-- §§ [/zzParity]

-- | X-rotation on qubit 0, identity on qubit 1.
--   Constructed by conjugating a 2-qubit diagonal by H⊗I.
xLow2 :: Natural -> Proc
xLow2 k =
  conj (H ⊗ I) ((Id ⊕ delta k) ⊕ (Id ⊕ delta k))

--------------------------------------------------------------------------------
-- Small examples

-- §§ [zxK]
zxK :: Natural -> Natural -> Proc
zxK kz kx = zK kz ∥ xK kx
-- §§ [/zxK]

testZX :: Natural -> IO ()
testZX n = runAndCompare True "test.qasm" n (zxK 2 3)

testZXPerms :: Natural -> IO ()
testZXPerms n = comparePerms "test.qasm" n (zxK 2 3)

-- §§ [zxyK]
zxyK :: Proc
zxyK = zK 1 ∥ xK 2 ∥ yK 3
-- §§ [/zxyK]

testZXYPerms :: Natural -> IO ()
testZXYPerms n = comparePerms "test.qasm" n zxyK

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
