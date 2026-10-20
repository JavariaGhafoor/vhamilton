-- | Energies of the chemistry ansatz written as ONE parallel VHamilton block,
-- compiled at a chosen Trotter level n through the ordinary pipeline
-- (Compiler + Qasm + the Qiskit statevector oracle).
--
--   G(α) = α1 YI ∥ α2 (IY + ZY)/2 ∥ α3 (IY − ZY)/2
--
-- At n = 0 the canonical schedule runs the three generators one after
-- another, which is exactly the sequential ansatz of ChemVQE; as n grows the
-- compiled circuit approaches the block's simultaneous evolution
-- e^{-2πi G(α)}.  YI anticommutes with ZY, so n matters.
--
-- Schedule "canon" is Compiler.compileTrotterWith.  Schedule "sym" alternates
-- the canonical schedule and its reverse between the 2^n repetitions: the
-- same atoms, the same count, a symmetric (second-order) product.
--
-- Ansatz "pair" writes one Pauli generator per clock: the noncommuting pair
-- [α1 YI ∥ α2 ZY] compiled at level n, then [α3 IY] compiled once at level 0
-- (IY commutes with both, so it needs no Trotterization).
--
-- Used to validate scripts/trotter_search/vh_numpy.py.
--
-- Usage:  trotter-chem-energies controlled|pair K R N canon|sym  < points
--   each stdin line: q1 q2 q3   (clock times q_j / 2^R)
-- Prints one line per point: q1 q2 q3 energy
module Main where

import Numeric.Natural (Natural)
import System.Environment (getArgs)
import Text.Printf (printf)

import ExamplesChem (fetchDecomp, chemProcWithIdentity, pauliParamProc, pauliSignedParamProc)
import Syntax (Proc, QCircuit, Valuation, (∥))
import Semantics (canonicalSchedule)
import Compiler (compileScheduleWith, compileTrotterWith)
import VQE (VAnsatz(..), VBlock(..), valuationForParams, energyForCircuit)

parallelBlock :: Proc
parallelBlock =
  pauliParamProc "YI" 0 "alpha1"
  ∥ pauliParamProc "IY" 1 "alpha2"
  ∥ pauliParamProc "ZY" 1 "alpha2"
  ∥ pauliParamProc "IY" 1 "alpha3"
  ∥ pauliSignedParamProc "ZY" 1 False "alpha3"

pairBlock, freeBlock :: Proc
pairBlock = pauliParamProc "YI" 0 "alpha1" ∥ pauliParamProc "ZY" 0 "alpha2"
freeBlock = pauliParamProc "IY" 0 "alpha3"

ansatzAt :: Natural -> VAnsatz
ansatzAt r = VAnsatz
  { vaParamNames = ["alpha1", "alpha2", "alpha3"]
  , vaPrecision  = r
  , vaBlocks     = [VBlock parallelBlock 0]
  }

-- | One block at level n, with the canonical or the alternating schedule.
compileBlock :: String -> Valuation -> Natural -> Proc -> QCircuit
compileBlock sched rho n p =
  case sched of
    "canon" -> compileTrotterWith rho n p
    "sym"   -> mconcat [ compileScheduleWith rho (if even i then s else reverse s)
                       | i <- [0 .. reps - 1] ]
    _       -> error "schedule must be canon or sym"
  where
    s    = canonicalSchedule p n
    reps = 2 ^ (fromIntegral n :: Int) :: Int

compileAt :: String -> String -> Natural -> Natural -> [Integer] -> QCircuit
compileAt ans sched r n qs =
  case ans of
    "controlled" -> compileBlock sched rho n parallelBlock
    "pair"       -> compileBlock sched rho n pairBlock <> compileTrotterWith rho 0 freeBlock
    _            -> error "ansatz must be controlled or pair"
  where
    rho = valuationForParams (ansatzAt r) qs

main :: IO ()
main = do
  args <- getArgs
  case args of
    [ans, k, r, n, sched] -> do
      d <- fetchDecomp "upper" (read k)
      let target = chemProcWithIdentity d
      pts <- map (map read . words) . filter (not . null) . lines <$> getContents
      mapM_ (\qs -> do
               (e, _) <- energyForCircuit target (compileAt ans sched (read r) (read n) qs)
               printf "%s %.15f\n" (unwords (map show qs)) e)
            pts
    _ -> error "usage: trotter-chem-energies controlled|pair K R N canon|sym < points"
