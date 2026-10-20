-- | Minimal noiseless VQE example implemented directly in VHamilton.
--
-- Target: H = X + Z, whose exact ground energy is -sqrt(2).
-- Ansatz family: exp(-2*pi*i * alpha * Y)|0>, represented by a single open
-- VHamilton process whose delta superscripts contain the free variable alpha.
-- The classical grid uses alpha = q/16, q = 0,...,7.
module Main where

import Text.Printf (printf)

import Syntax
  ( Clifford(..), (·)
  , Proc(..), (∥)
  , BareAction(..), deltaPow, (⊕), Action(..)
  )
import VQE
  ( VBlock(..), VAnsatz(..), VQEResult(..), gridSearch )

act :: BareAction -> Proc
act = Atom . Bare

conj :: Clifford -> BareAction -> Proc
conj u d = Atom (Conj u (Bare d))

-- Exact closed Pauli generators for the target Hamiltonian.
pauliZ :: Proc
pauliZ = act (deltaPow 0 1 ⊕ deltaPow 0 (-1))

pauliX :: Proc
pauliX = conj H (deltaPow 0 1 ⊕ deltaPow 0 (-1))

-- Open VHamilton family alpha * Y.
pauliYVar :: Proc
pauliYVar =
  conj (S · H)
    (DeltaParam 0 "alpha" ⊕ NegDeltaParam 0 "alpha")

main :: IO ()
main = do
  let target = pauliX ∥ pauliZ
      ansatz = VAnsatz
        { vaParamNames = ["alpha"]
        , vaPrecision = 4            -- rho(alpha) = q / 16
        , vaBlocks = [VBlock pauliYVar 0]
        }
      grids = [[0..7]]

  (best, allResults) <- gridSearch target ansatz grids

  putStrLn "VHamilton VQE: H = X + Z"
  putStrLn "Open process: Y(alpha), with rho(alpha)=q/16 and q in {0,...,7}"
  putStrLn "Ansatz: exp(-2*pi*i*alpha*Y)|0>"
  putStrLn ""
  mapM_ (\r ->
    putStrLn $ "q=" ++ show (head (vrParams r))
            ++ "  E=" ++ printf "%.12f" (vrEnergy r))
    allResults

  putStrLn ""
  putStrLn $ "best q       = " ++ show (vrParams best)
  putStrLn $ "best energy  = " ++ printf "%.12f" (vrEnergy best)
  putStrLn $ "exact ground = " ++ printf "%.12f" (vrGroundEnergy best)
