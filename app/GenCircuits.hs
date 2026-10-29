-- | Generate quantikz circuit diagrams.
--
-- Run from the code/ directory:
--   cabal run gen-circuits -- <output-dir>
module Main where

import System.Environment (getArgs)
import System.Directory (createDirectoryIfMissing)
import System.FilePath ((</>))

import Compiler (compileTrotter)
import Quantikz (renderQuantikz)
import Examples (zxK, zzParity)

-- | Wrap a quantikz environment in a minimal standalone document so it can
-- be compiled to a self-contained PDF and included via \includegraphics,
-- keeping quantikz out of the main document (avoids tikz/math font conflicts).
standalone :: String -> String
standalone body = unlines
  [ "\\documentclass{standalone}"
  , "\\usepackage{quantikz}"
  , "\\begin{document}"
  , body
  , "\\end{document}"
  ]

writeCircuit :: FilePath -> String -> IO ()
writeCircuit path tex =
  writeFile path (standalone tex) >> putStrLn ("  wrote " ++ path)

main :: IO ()
main = do
  args <- getArgs
  let (outDir, which) = case args of
        [d]    -> (d, "all")
        [d, w] -> (d, w)
        _      -> ("generated", "all")
  createDirectoryIfMissing True outDir

  let genZX = do
        -- ZX at refinement level 0: one Z-step then one X-step
        writeCircuit (outDir </> "circuit_zx_n0.tex") $
          renderQuantikz (compileTrotter 0 (zxK 2 3))
        -- ZX at refinement level 1: two repetitions of half-angle steps
        writeCircuit (outDir </> "circuit_zx_n1.tex") $
          renderQuantikz (compileTrotter 1 (zxK 2 3))
      genZZ =
        writeCircuit (outDir </> "circuit_zzparity_n0.tex") $
          renderQuantikz (compileTrotter 0 (zzParity 4))

  putStrLn "Generating circuits..."
  case which of
    "all" -> genZX >> genZZ >> putStrLn "All circuits written."
    "zx"  -> genZX >> putStrLn "ZX circuits written."
    "zz"  -> genZZ >> putStrLn "ZZ circuit written."
    _ -> error $ "unknown circuit experiment '" ++ which ++ "' (use all, zx, or zz)"
