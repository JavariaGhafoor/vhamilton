-- | One-layer QAOA for MaxCut on a single edge, written in VHamilton.
--
-- Target (cost) Hamiltonian:  C = Z⊗Z, ground energy -1 (the two spins
-- disagree, i.e. the edge is cut).
--
-- Ansatz, as three sequential blocks (sequencing lives in the VQE layer;
-- ∥ inside a block is still Hamiltonian addition):
--
--   prep      = YI/8 ∥ IY/8          closed Hamilton block, |00> ↦ |++>
--   cost(γ)   = γ ZZ                 one open clock on the coupling
--   mixer(β)  = β XI ∥ β IX          one shared clock on both fields
--
-- At search precision r the valuation is rho(α) = q/2^r, so 'DeltaParam 0 α'
-- under rho is the paper's δ_r^α with ν(α) = q.  Since e^{-2πiγZZ} and
-- e^{-2πiβ(XI+IX)} are periodic (up to global phase) with period 1/2 in
-- γ and β, the grid q ∈ {0,…,2^(r-1)-1} covers one full period.
--
-- Usage:  qaoa-demo R1 [R2 ...]   (default: 3 4 5)
-- Writes  qaoa_grid_r<R>.csv  with columns q_gamma,q_beta,energy.
module Main where

import Numeric.Natural (Natural)
import System.Environment (getArgs)
import Text.Printf (printf)

import Syntax (Proc, (∥))
import ExamplesChem (ChemTerm(..), chemTerm, pauliProc, pauliParamProc)
import VQE (VBlock(..), VAnsatz(..), VQEResult(..), gridSearch)

-- | Closed Hamilton block (YI + IY)/8: rotates |00> to |++>.
prepPlus :: Proc
prepPlus = chemTerm (ChemTerm "YI" 1 3) ∥ chemTerm (ChemTerm "IY" 1 3)

-- | Open cost block γ ZZ.
costBlock :: Proc
costBlock = pauliParamProc "ZZ" 0 "gamma"

-- | Open mixer block β (XI + IX): one clock shared by two interactions.
mixerBlock :: Proc
mixerBlock = pauliParamProc "XI" 0 "beta" ∥ pauliParamProc "IX" 0 "beta"

ansatzAt :: Natural -> VAnsatz
ansatzAt r = VAnsatz
  { vaParamNames = ["gamma", "beta"]
  , vaPrecision  = r
  , vaBlocks     = [VBlock prepPlus 0, VBlock costBlock 0, VBlock mixerBlock 0]
  }

runAt :: Natural -> IO ()
runAt r = do
  let half = 2 ^ (fromIntegral r - 1 :: Int) - 1 :: Integer
      grid = [[0 .. half], [0 .. half]]
  (best, results) <- gridSearch (pauliProc "ZZ") (ansatzAt r) grid
  let csv = "qaoa_grid_r" ++ show r ++ ".csv"
  writeFile csv $ unlines $
    "q_gamma,q_beta,energy" :
    [ printf "%d,%d,%.12f" g b (vrEnergy res)
    | res <- results, (g, b) <- pair (vrParams res) ]
  printf "r=%d  grid=%d  best q=%s  best E=%.12f  exact ground=%.12f  -> %s\n"
    (fromIntegral r :: Int) (length results) (show (vrParams best))
    (vrEnergy best) (vrGroundEnergy best) csv

pair :: [Integer] -> [(Integer, Integer)]
pair [g, b] = [(g, b)]
pair qs     = error ("qaoa-demo: expected two parameters, got " ++ show qs)

main :: IO ()
main = do
  args <- getArgs
  let rs = if null args then [3, 4, 5] else map read args
  mapM_ runAt rs
