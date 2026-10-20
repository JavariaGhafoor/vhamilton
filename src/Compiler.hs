-- | Compilation pipeline from processes to circuits (§Implementation, fig:circuits).
--
--   Proc  --(canonicalSchedule)-->  Schedule  --(this module)-->  QCircuit
--
-- compileClifford  -- Clifford u to gate sequence (offset-aware, structural)
-- diagonalAngles   -- phase angles of ⟦d⟧_U = exp(-2πi ⟦d⟧_H) for bare actions (fig:sem)
-- compileAction    -- action a to QCircuit, applying conjugation for Conj (fig:sem)
-- compileTrotter   -- Trotter approximation: (⟦sch(p)/2^n⟧_U)^{2^n} → ⟦p⟧_U (eq:trotter)
module Compiler where

import Numeric.Natural (Natural)
import qualified Data.Map.Strict as M

-- Open VHamilton terms are compiled under an explicit valuation.
import Syntax
  ( Clifford(..), dimClifford, daggerC
  , Proc, Valuation, parameterValue
  , BareAction(..), Action(..), Schedule, dimBareAction
  , QGate(..), QCircuit(..), oneC
  )
import Semantics (canonicalSchedule)

-- | Convert a matrix dimension (power of 2) to a qubit count for gate addressing.
ilog2 :: Natural -> Natural
ilog2 1 = 0
ilog2 n = 1 + ilog2 (n `div` 2)

--------------------------------------------------------------------------------
-- Clifford -> QCircuit

compileClifford :: Natural -> Clifford -> QCircuit
compileClifford off = \case
  I       -> mempty
  H       -> oneC (QH off)
  S       -> oneC (QS off)
  CX      -> oneC (QCX  off (off + 1))
  SWAP    -> oneC (QSWAP off (off + 1))
  u :⊗: v ->
    compileClifford off u <>
    compileClifford (off + ilog2 (dimClifford u)) v
  -- Matrix semantics is [[u · v]] = [[u]][[v]].  Since circuit gates
  -- execute chronologically left-to-right, emit v first and then u so the
  -- resulting circuit matrix is U V.
  u :⋅: v ->
    compileClifford off v <>
    compileClifford off u

--------------------------------------------------------------------------------
-- Diagonal phase angles
--
-- A bare action's Hamiltonian value is converted to a phase -2π times
-- that value.  The precision index k is still halved by schedule refinement.

phaseAngle :: Natural -> Double
phaseAngle k =
  -(2 * pi) /
  fromIntegral (2 ^ (fromIntegral k :: Int) :: Integer)

-- §§ [diag-angles]
diagonalAnglesWith :: Valuation -> BareAction -> [Double]
diagonalAnglesWith rho = \case
  Id                 -> [0]
  Delta k            -> [phaseAngle k]
  NegDelta k         -> [-phaseAngle k]
  DeltaPow k j       -> [fromInteger j * phaseAngle k]
  DeltaParam k a     -> [parameterValue rho a * phaseAngle k]
  NegDeltaParam k a  -> [-parameterValue rho a * phaseAngle k]
  a :⊕: b            -> diagonalAnglesWith rho a ++ diagonalAnglesWith rho b

-- Closed Hamilton programs need no valuation; open VHamilton programs do.
diagonalAngles :: BareAction -> [Double]
diagonalAngles = diagonalAnglesWith M.empty
-- §§ [/diag-angles]

--------------------------------------------------------------------------------
-- Action -> QCircuit

compileActionWith :: Valuation -> Natural -> Action -> QCircuit
compileActionWith rho off = \case
  Bare d ->
    emitDiagonal off (ilog2 (dimBareAction d)) (diagonalAnglesWith rho d)
  -- Emit U†, then the action, then U, in chronological order.
  Conj u a ->
    compileClifford off (daggerC u) <>
    compileActionWith rho off a <>
    compileClifford off u

compileAction :: Natural -> Action -> QCircuit
compileAction = compileActionWith M.empty

emitDiagonal :: Natural -> Natural -> [Double] -> QCircuit
emitDiagonal off n phis = case (n, phis) of
  (0, [phi]) -> oneC (QGPhase phi)
  (0, _)     ->
    error "emitDiagonal: 0-qubit diagonal must have exactly one phase"
  _          -> oneC (QDiagonalAt off n phis)

--------------------------------------------------------------------------------
-- Schedule -> QCircuit

-- | Compile a schedule with or without an explicit symbolic valuation.
compileScheduleWith :: Valuation -> Schedule -> QCircuit
compileScheduleWith rho = foldMap (compileActionWith rho 0)

compileSchedule :: Schedule -> QCircuit
compileSchedule = compileScheduleWith M.empty

-- §§ [compile-trotter]
compileTrotterWith :: Valuation -> Natural -> Proc -> QCircuit
compileTrotterWith rho n proc =
  mconcat (replicate reps step)
  where
    step = compileScheduleWith rho (canonicalSchedule proc n)
    reps = 2 ^ (fromIntegral n :: Int)

compileTrotter :: Natural -> Proc -> QCircuit
compileTrotter = compileTrotterWith M.empty
-- §§ [/compile-trotter]

-- | Compile an arbitrary schedule repeated 2^n times.
compileScheduleN :: Natural -> Schedule -> QCircuit
compileScheduleN n sched =
  mconcat (replicate reps (compileSchedule sched))
  where
    reps = 2 ^ (fromIntegral n :: Int)

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
