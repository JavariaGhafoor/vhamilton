-- | Syntax and dimensional type system for Hamilton and minimal VHamilton.
--
--   u ::= I | H | S | CX | SWAP | u · u | u ⊗ u
--   d ::= id | δ_k | -δ_k | δ_k^j | δ_k^alpha_j | -δ_k^alpha_j | d ⊕ d
--   a ::= d | u a u†
--   p ::= a | p ∥ p
--
-- The only VHamilton additions to Hamilton are the two parameterized bare
-- actions. No general scalar expressions, scalar type, or scaling constructor.
-- A valuation rho maps a parameter name to its real (usually dyadic) time.
module Syntax where

import Numeric.Natural (Natural)
import Data.Text (Text, pack)
import qualified Data.Text as Text
import GHC.Generics (Generic)
import qualified Data.Aeson as A
import qualified Data.Map.Strict as M
import qualified Data.Set as SSet

--------------------------------------------------------------------------------
-- Valuations live outside the object language; they are supplied by VQE.

type Valuation = M.Map String Double

-- | Look up a symbolic parameter without changing the process syntax.
parameterValue :: Valuation -> String -> Double
parameterValue rho a =
  case M.lookup a rho of
    Just x  -> x
    Nothing -> error ("unbound VHamilton parameter: " ++ a)

--------------------------------------------------------------------------------
-- Clifford transformations
--   u  ::=  I | H | S | CX | SWAP | u ⊗ u | u · u

infixr 7 :⊗:
infixr 6 :⋅:

data Clifford
  = I | H | S | CX | SWAP
  | Clifford :⋅: Clifford
  | Clifford :⊗: Clifford
  deriving (Eq, Ord, Show, Generic)

infixr 7 ⊗
infixr 6 ·

(⊗) :: Clifford -> Clifford -> Clifford
(⊗) = (:⊗:)

(·) :: Clifford -> Clifford -> Clifford
(·) = (:⋅:)

-- | Adjoint of a Clifford.
daggerC :: Clifford -> Clifford
daggerC = \case
  I       -> I
  H       -> H
  S       -> S · S · S
  CX      -> CX
  SWAP    -> SWAP
  u :⊗: v -> daggerC u ⊗ daggerC v
  u :⋅: v -> daggerC v · daggerC u

-- | Matrix dimension of a Clifford (paper types: I:2, H:2, S:2, CX:4,
-- SWAP:4; u⊗v:nm; u·v:n with n=m enforced).
dimClifford :: Clifford -> Natural
dimClifford = \case
  I       -> 2
  H       -> 2
  S       -> 2
  CX      -> 4
  SWAP    -> 4
  u :⊗: v -> dimClifford u * dimClifford v
  u :⋅: v -> arityEq "Clifford (·)" (dimClifford u) (dimClifford v)

prettyClifford :: Clifford -> String
prettyClifford = \case
  I       -> "I"
  H       -> "H"
  S       -> "S"
  CX      -> "CX"
  SWAP    -> "SWAP"
  u :⊗: v -> prettyClifford u ++ " ⊗ " ++ prettyClifford v
  u :⋅: v -> prettyClifford u ++ " · " ++ prettyClifford v

instance A.ToJSON Clifford where
  toJSON I         = tag "I"
  toJSON H         = tag "H"
  toJSON S         = tag "S"
  toJSON CX        = tag "CX"
  toJSON SWAP      = tag "SWAP"
  toJSON (u :⊗: v) = tagC "Tensor"  [A.toJSON u, A.toJSON v]
  toJSON (u :⋅: v) = tagC "Compose" [A.toJSON u, A.toJSON v]

--------------------------------------------------------------------------------
-- Bare actions: original Hamilton + exactly two parameterized constructors.
-- [[δ_k]] = 2^-k; [[-δ_k]] = -2^-k; [[δ_k^j]] = j * 2^-k.
-- [[±δ_k^alpha]]_rho = ±rho(alpha) * 2^-k.

infixl 6 :⊕:

data BareAction
  = Id
  | Delta Natural
  | NegDelta Natural
  | DeltaPow Natural Integer
  | DeltaParam Natural String
  | NegDeltaParam Natural String
  | BareAction :⊕: BareAction
  deriving (Eq, Ord, Show)

infixl 6 ⊕

(⊕) :: BareAction -> BareAction -> BareAction
(⊕) = (:⊕:)

-- Familiar Haskell helpers, each returning a genuine bare action.
delta :: Natural -> BareAction
delta = Delta

negDelta :: Natural -> BareAction
negDelta = NegDelta

deltaPow :: Natural -> Integer -> BareAction
deltaPow = DeltaPow

deltaParam :: Natural -> String -> BareAction
deltaParam = DeltaParam

negDeltaParam :: Natural -> String -> BareAction
negDeltaParam = NegDeltaParam

dimBareAction :: BareAction -> Natural
dimBareAction = \case
  Id               -> 1
  Delta _          -> 1
  NegDelta _       -> 1
  DeltaPow _ _     -> 1
  DeltaParam _ _   -> 1
  NegDeltaParam _ _ -> 1
  a :⊕: b          -> dimBareAction a + dimBareAction b

prettyBareAction :: BareAction -> String
prettyBareAction = \case
  Id                 -> "id"
  Delta k            -> "δ" ++ show k
  NegDelta k         -> "-δ" ++ show k
  DeltaPow k j       -> "δ" ++ show k ++ "^" ++ show j
  DeltaParam k a     -> "δ" ++ show k ++ "^{" ++ a ++ "}"
  NegDeltaParam k a  -> "-δ" ++ show k ++ "^{" ++ a ++ "}"
  a :⊕: b            -> "(" ++ prettyBareAction a ++ " ⊕ " ++ prettyBareAction b ++ ")"

instance A.ToJSON BareAction where
  toJSON Id                   = tag "Id"
  toJSON (Delta k)            = tagV "Delta" (toInteger k)
  toJSON (NegDelta k)         = tagV "NegDelta" (toInteger k)
  toJSON (DeltaPow k j)       = tagC "DeltaPow" [A.toJSON (toInteger k), A.toJSON j]
  toJSON (DeltaParam k a)     = tagC "DeltaParam" [A.toJSON (toInteger k), A.toJSON a]
  toJSON (NegDeltaParam k a)  = tagC "NegDeltaParam" [A.toJSON (toInteger k), A.toJSON a]
  toJSON (a :⊕: b)            = tagC "Oplus" [A.toJSON a, A.toJSON b]

--------------------------------------------------------------------------------
-- Actions
--   a  ::=  d | u a u†

data Action
  = Bare BareAction
  | Conj Clifford Action
  deriving (Eq, Ord, Show)

-- | Matrix dimension of an action; enforces u:n and a:n agree in u a u†.
dimAction :: Action -> Natural
dimAction = \case
  Bare d    -> dimBareAction d
  Conj u a  -> arityEq "Conj (u a u†)" (dimClifford u) (dimAction a)

prettyAction :: Action -> String
prettyAction = \case
  Bare d    -> prettyBareAction d
  Conj u a  -> "Conj(" ++ prettyClifford u ++ ", " ++ prettyAction a ++ ")"

instance A.ToJSON Action where
  toJSON (Bare d)   = tagV "Base" d
  toJSON (Conj u a) = tagC "Conj" [A.toJSON u, A.toJSON a]

--------------------------------------------------------------------------------
-- Processes
--   p  ::=  a | p ∥ p

infixr 4 :∥:

-- §§ [proc-data]
data Proc
  = Atom Action
  | Proc :∥: Proc
-- §§ [/proc-data]
  deriving (Eq, Ord, Show, Generic)

infixr 4 ∥

(∥) :: Proc -> Proc -> Proc
(∥) = (:∥:)

-- | Matrix dimension of a process; enforces both sides of ∥ act on the same space.
dimProc :: Proc -> Natural
dimProc = \case
  Atom a     -> dimAction a
  p :∥: q    -> arityEq "Proc (∥)" (dimProc p) (dimProc q)

prettyProc :: Proc -> String
prettyProc = go (0 :: Int)
  where
    go _ (Atom a)  = prettyAction a
    go p (q :∥: r) =
      paren (p > 4) (go 4 q ++ " ∥ " ++ go 4 r)
    paren True  s  = "(" ++ s ++ ")"
    paren False s  = s

instance A.ToJSON Proc where
  toJSON (Atom a)  = tagV "Base" a
  toJSON (p :∥: q) = tagC "Par" [A.toJSON p, A.toJSON q]

--------------------------------------------------------------------------------
-- VHamilton parameter utilities

-- | Names occurring in the two (and only two) parameterized bare actions.
freeParams :: Proc -> [String]
freeParams = SSet.toAscList . goP
  where
    goP (Atom a)   = goA a
    goP (p :∥: q)  = goP p `SSet.union` goP q

    goA (Bare d)   = goD d
    goA (Conj _ a) = goA a

    goD Id                    = SSet.empty
    goD (Delta _)             = SSet.empty
    goD (NegDelta _)          = SSet.empty
    goD (DeltaPow _ _)        = SSet.empty
    goD (DeltaParam _ a)      = SSet.singleton a
    goD (NegDeltaParam _ a)   = SSet.singleton a
    goD (a :⊕: b)             = goD a `SSet.union` goD b

--------------------------------------------------------------------------------
-- Schedules: lists of actions

type Schedule = [Action]

prettySchedule :: Schedule -> Text
prettySchedule acts = Text.concat
  [ Text.show i `Text.append` ": " `Text.append` (pack $ prettyAction a)
  | (i, a) <- zip [(0 :: Int)..] acts
  ]

--------------------------------------------------------------------------------
-- QASM gate types and circuits

data QGate
  = QGPhase Double
  | QH    Natural
  | QS    Natural
  | QX    Natural
  | QCX   Natural Natural
  | QSWAP Natural Natural
  | QP    Double Natural
  | QCP   Double Natural Natural
  | QDiagonalAt Natural Natural [Double]
  deriving (Eq, Show)

gateQubits :: QGate -> [Natural]
gateQubits = \case
  QGPhase _           -> []
  QH i                -> [i]
  QS i                -> [i]
  QX i                -> [i]
  QCX c t             -> [c, t]
  QSWAP a b           -> [a, b]
  QP _ i              -> [i]
  QCP _ c t           -> [c, t]
  QDiagonalAt off n _ ->
    if n == 0 then [] else [off .. off + n - 1]

newtype QCircuit = QCircuit { getQCircuit :: [QGate] }
  deriving (Eq, Show)

instance Semigroup QCircuit where
  QCircuit a <> QCircuit b = QCircuit (a ++ b)

instance Monoid QCircuit where
  mempty  = QCircuit []
  mappend = (<>)

oneC :: QGate -> QCircuit
oneC g = QCircuit [g]

-- | Minimum qubit register size for a circuit.
numQubitsQC :: QCircuit -> Natural
numQubitsQC (QCircuit gs) =
  case concatMap gateQubits gs of
    [] -> 0
    qs -> maximum qs + 1

--------------------------------------------------------------------------------
-- JSON helpers (used by ToJSON instances above)

tag :: String -> A.Value
tag t = A.object ["tag" A..= t]

tagV :: A.ToJSON v => String -> v -> A.Value
tagV t v = A.object ["tag" A..= t, "contents" A..= v]

tagC :: String -> [A.Value] -> A.Value
tagC t vs = A.object ["tag" A..= t, "contents" A..= vs]

--------------------------------------------------------------------------------
-- Internal helper

arityEq :: String -> Natural -> Natural -> Natural
arityEq ctx a b
  | a == b    = a
  | otherwise = error
      (ctx ++ ": dimension mismatch (" ++ show a ++ " vs " ++ show b ++ ")")
