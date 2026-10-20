-- | VQE over a noncommuting VHamilton block, compiled at several Trotter levels.
--
-- Target:  H = X + Y + Z  (one qubit, ground energy -sqrt 3).
-- Ansatz:  one open block  Z(alpha) ∥ X(beta)  acting on |0>.
--
-- The block's two terms do not commute, so its unit-time evolution
-- e^{-2πi(αZ+βX)} is compiled by Lie-Trotter interleaving at level n.
-- Each level n gives a *different* circuit family; only as n grows do the
-- compiled energies converge to those of the program's own trajectories.
-- The experiment asks what an optimizer finds when it only ever sees the
-- compiled energies at a fixed n.
--
-- At precision r the valuation is rho(α) = q/2^r with q ∈ {0,…,2^r-1}.
--
-- Usage:  trotter-vqe-demo R N1 [N2 ...]
-- Writes  trotter_grid_r<R>_n<N>.csv with columns q_alpha,q_beta,energy.
module Main where

import Numeric.Natural (Natural)
import System.Environment (getArgs)
import Text.Printf (printf)

import Syntax
  ( Clifford(..), (·), Proc(..), (∥)
  , BareAction(..), deltaPow, (⊕), Action(..) )
import VQE (VBlock(..), VAnsatz(..), VQEResult(..), gridSearch)

act :: BareAction -> Proc
act = Atom . Bare

conj :: Clifford -> BareAction -> Proc
conj u d = Atom (Conj u (Bare d))

zFixed, xFixed, yFixed :: Proc
zFixed = act (deltaPow 0 1 ⊕ deltaPow 0 (-1))
xFixed = conj H (deltaPow 0 1 ⊕ deltaPow 0 (-1))
yFixed = conj (S · H) (deltaPow 0 1 ⊕ deltaPow 0 (-1))

-- | Open block  alpha Z ∥ beta X.
block :: Proc
block =
  act (DeltaParam 0 "alpha" ⊕ NegDeltaParam 0 "alpha")
    ∥ conj H (DeltaParam 0 "beta" ⊕ NegDeltaParam 0 "beta")

runAt :: Natural -> Natural -> IO ()
runAt r n = do
  let top = 2 ^ (fromIntegral r :: Int) - 1 :: Integer
      ansatz = VAnsatz
        { vaParamNames = ["alpha", "beta"]
        , vaPrecision  = r
        , vaBlocks     = [VBlock block n]
        }
  (best, results) <- gridSearch (xFixed ∥ yFixed ∥ zFixed) ansatz [[0 .. top], [0 .. top]]
  let csv = printf "trotter_grid_r%d_n%d.csv" (fromIntegral r :: Int) (fromIntegral n :: Int) :: String
  writeFile csv $ unlines $
    "q_alpha,q_beta,energy" :
    [ printf "%d,%d,%.12f" a b (vrEnergy res)
    | res <- results, (a, b) <- pair (vrParams res) ]
  printf "r=%d n=%d  grid=%d  best q=%s  best E=%.12f  exact ground=%.12f  -> %s\n"
    (fromIntegral r :: Int) (fromIntegral n :: Int) (length results)
    (show (vrParams best)) (vrEnergy best) (vrGroundEnergy best) csv

pair :: [Integer] -> [(Integer, Integer)]
pair [a, b] = [(a, b)]
pair qs     = error ("trotter-vqe-demo: expected two parameters, got " ++ show qs)

main :: IO ()
main = do
  args <- getArgs
  case map read args of
    (r : ns@(_ : _)) -> mapM_ (runAt r) ns
    _ -> error "usage: trotter-vqe-demo R N1 [N2 ...]"
