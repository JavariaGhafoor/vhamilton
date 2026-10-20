-- | Three minimal VHamilton implementation examples.
--
-- 1. One open process, two valuations.
-- 2. A non-commuting variational process with ordinary Trotter convergence.
-- 3. A commuting variational ZZ interaction, exact at every refinement.
module Main where

import qualified Data.Map.Strict as M
import Text.Printf (printf)

import Syntax (freeParams, prettyProc)
import ExamplesVHamilton (zV, zxV, zzV, instantiateAndCompare)

main :: IO ()
main = do
  putStrLn "VHamilton minimal examples"
  putStrLn ""

  let p1 = zV "alpha"
  putStrLn "1. A symbolic process evaluated at two different times"
  putStrLn $ "   p = " ++ prettyProc p1
  putStrLn $ "   free parameters = " ++ show (freeParams p1)
  mapM_ (showValuation p1) [1/8, 3/8]
  putStrLn ""

  putStrLn "2. Non-commuting family: alpha Z || beta X"
  let rhoZX = M.fromList [("alpha", 1/4), ("beta", 1/8)]
  mapM_ (showFrob rhoZX zxV) [0..4]
  putStrLn ""

  putStrLn "3. Commuting/single-interaction family: gamma ZZ"
  let rhoZZ = M.fromList [("gamma", 3/16)]
  mapM_ (showFrob rhoZZ zzV) [0..3]

showValuation p a = do
  let rho = M.fromList [("alpha", a)]
  e <- instantiateAndCompare rho 0 p
  putStrLn $ "   rho(alpha)=" ++ printf "%.4f" a
          ++ " in " ++ prettyProc p
          ++ ", Frobenius=" ++ printf "%.3e" e

showFrob rho p n = do
  e <- instantiateAndCompare rho n p
  putStrLn $ "   refinement n=" ++ show n
          ++ "  Frobenius=" ++ printf "%.6f" e
