-- | Generate LaTeX table fragments (§Examples tables and §Application tables).
--
-- Run from the code/ directory:
--   cabal run gen-tables -- <output-dir>
-- Writes one .tex file per table into <output-dir>.
--
-- genZXConv      -- ZX Trotter convergence table (§Examples)
-- genZXOrd       -- ZX schedule orderings (§Examples, Rule R1)
-- genZXYOrd      -- ZXY schedule orderings (§Examples, Rule R1)
-- genZZParity    -- exact result for commuting ZZ (§Examples)
-- genChemTrotter -- DMANH+ Trotter convergence sweep (§Application)
-- genDyadicError -- dyadic approximation error vs. level k (§Application)
module Main where

import Numeric.Natural (Natural)
import System.Environment (getArgs)
import System.Directory (createDirectoryIfMissing)
import System.FilePath ((</>))
import Control.Monad (forM)
import Data.List (sortOn, intercalate)
import Text.Printf (printf)

import Syntax (Proc, Schedule, Action(..), Clifford(..), prettyClifford)
import Compiler (compileTrotter, compileScheduleN)
import Semantics (schedulePerms)
import Denot (compareToIdeal, compareIdeals)
import Examples (zxK, zxyK, zzParity)
import ExamplesChem (fetchDecomp, chemProc)

--------------------------------------------------------------------------------

frobTrotter :: Natural -> Proc -> IO Double
frobTrotter n p = fst <$> compareToIdeal p (compileTrotter n p)

frobPerm :: Natural -> Proc -> Schedule -> IO Double
frobPerm n p sched = fst <$> compareToIdeal p (compileScheduleN n sched)

writeTable :: FilePath -> String -> [String] -> [String] -> IO ()
writeTable path colSpec header rows =
  writeFile path $ unlines $
    [ "\\begin{tabular}{" ++ colSpec ++ "}"
    , "\\toprule"
    , intercalate " & " header ++ " \\\\"
    , "\\midrule"
    ] ++
    map (\r -> r ++ " \\\\") rows ++
    [ "\\bottomrule"
    , "\\end{tabular}"
    ]

--------------------------------------------------------------------------------
-- Action labels (derived from the Clifford wrapper on each action)

labelAction :: Action -> String
labelAction (Bare _)           = "Z"
labelAction (Conj H _)         = "X"
labelAction (Conj (S :⋅: H) _) = "Y"
labelAction (Conj u _)         = prettyClifford u

labelSchedule :: Schedule -> String
labelSchedule = intercalate "→" . map labelAction

--------------------------------------------------------------------------------
-- ZX convergence (n = 0..4)

genZXConv :: FilePath -> IO ()
genZXConv outDir = do
  putStrLn "Generating ZX convergence table..."
  let p = zxK 2 3
  rows <- forM [0..4] $ \n -> do
    frob <- frobTrotter n p
    return $ show n ++ " & " ++ show (2^n :: Int) ++ " & " ++ printf "%.4f" frob
  writeTable (outDir </> "table_zx_conv.tex") "r r r"
    ["$n$", "reps $2^n$", "Frobenius"] rows
  putStrLn "  done."

--------------------------------------------------------------------------------
-- ZX orderings at n=1

genZXOrd :: FilePath -> IO ()
genZXOrd outDir = do
  putStrLn "Generating ZX orderings table..."
  let p     = zxK 2 3
      perms = schedulePerms p 1
  frobs <- mapM (frobPerm 1 p) perms
  let ranked = sortOn fst (zip frobs perms)
      rows   = map (\(frob, sched) ->
                  labelSchedule sched ++ " & " ++ printf "%.4f" frob)
               ranked
  writeTable (outDir </> "table_zx_ord.tex") "l r"
    ["ordering", "Frobenius"] rows
  putStrLn "  done."

--------------------------------------------------------------------------------
-- ZXY orderings at n=1

genZXYOrd :: FilePath -> IO ()
genZXYOrd outDir = do
  putStrLn "Generating ZXY orderings table..."
  let p     = zxyK
      perms = schedulePerms p 1
  frobs <- mapM (frobPerm 1 p) perms
  let ranked = sortOn fst (zip frobs perms)
      rows   = map (\(frob, sched) ->
                  labelSchedule sched ++ " & " ++ printf "%.4f" frob)
               ranked
  writeTable (outDir </> "table_zxy_ord.tex") "l r"
    ["ordering", "Frobenius"] rows
  putStrLn "  done."

--------------------------------------------------------------------------------
-- zzParity exact result

genZZParity :: FilePath -> IO ()
genZZParity outDir = do
  putStrLn "Generating zzParity exact result..."
  frob <- frobTrotter 0 (zzParity 4)
  writeFile (outDir </> "frob_zzparity.tex") (printf "%.4f" frob)
  putStrLn "  done."

--------------------------------------------------------------------------------
-- Chemistry Trotter sweep (level=12, n=0..8)

genChemTrotter :: FilePath -> IO ()
genChemTrotter outDir = do
  putStrLn "Generating chemistry Trotter table (this may take a while)..."
  d <- fetchDecomp "upper" 12
  let p = chemProc d
  rows <- forM [0..8] $ \n -> do
    frob <- frobTrotter n p
    putStrLn $ "  n=" ++ show n ++ " -> " ++ printf "%.4f" frob
    return $ show n ++ " & " ++ show (2^n :: Int) ++ " & " ++ printf "%.4f" frob
  writeTable (outDir </> "table_chem_trotter.tex") "r r r"
    ["$n$", "reps $2^n$", "Frobenius"] rows
  putStrLn "  done."

--------------------------------------------------------------------------------
-- Dyadic approximation error (compare level k to reference level 20)
-- Measures ‖exp(-2πi H_k) - exp(-2πi H_ref)‖_F as a function of k.

genDyadicError :: FilePath -> IO ()
genDyadicError outDir = do
  putStrLn "Generating dyadic error table..."
  dRef <- fetchDecomp "upper" 20   -- reference: 2^(-20) ≈ 10^(-6) vs errors ~ 0.001
  let pRef = chemProc dRef
  rows <- forM [3, 4, 6, 8, 10, 12] $ \k -> do
    dk <- fetchDecomp "upper" k
    let pk = chemProc dk
    frob <- compareIdeals pk pRef
    putStrLn $ "  k=" ++ show k ++ " -> " ++ printf "%.4f" frob
    return $ show k ++ " & " ++ printf "%.4f" frob
  writeTable (outDir </> "table_dyadic.tex") "r r"
    ["$k$", "dyadic error"] rows
  putStrLn "  done."

--------------------------------------------------------------------------------

main :: IO ()
main = do
  args <- getArgs
  let (outDir, which) = case args of
        [d]    -> (d, "all")
        [d, w] -> (d, w)
        _      -> ("generated", "all")
  createDirectoryIfMissing True outDir
  case which of
    "all"           -> do
      genZXConv      outDir
      genZXOrd       outDir
      genZXYOrd      outDir
      genZZParity    outDir
      genChemTrotter outDir
      genDyadicError outDir
      putStrLn "All tables written."
    "zx-conv"       -> genZXConv outDir
    "zx-orderings"  -> genZXOrd outDir
    "zxy-orderings" -> genZXYOrd outDir
    "zz-parity"     -> genZZParity outDir
    "chem-trotter"  -> genChemTrotter outDir
    "dyadic-error"  -> genDyadicError outDir
    _ -> error $ "unknown table experiment '" ++ which ++
                 "' (use all, zx-conv, zx-orderings, zxy-orderings, zz-parity, chem-trotter, or dyadic-error)"
