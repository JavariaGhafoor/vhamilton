-- | Halving, schedule extraction, and schedule halving.
--
-- Implements halving, schedule extraction and schedule halving, and Rule R1:
--   seg0 / seg1 / segProc  -- n-way successive halving p / 2^n
--   sch                    -- schedule of a process: sch(p)
--   segSchedule            -- schedule halving: [a] / 2^n
--   canonicalSchedule      -- sch(p) / 2^n, extraction and halving composed
--   schedulePerms          -- all permutations of sch(p) / 2^n         (Rule R1)
--
-- The canonical schedule repeated 2^n times is the Trotter circuit compiled
-- by Compiler.compileTrotter.
module Semantics where

import Numeric.Natural (Natural)
import Data.List (permutations)

import Syntax (Proc(..), BareAction(..), Action(..), Schedule)

--------------------------------------------------------------------------------
-- n-way successive halving: p / 2^n
--
--   id / 2^n          = id
--   (±)δ_k^(alpha) / 2^n = (±)δ_{k+n}^(alpha), and likewise for Hamilton generators
--   (d ⊕ d') / 2^n    = (d / 2^n) ⊕ (d' / 2^n)
--   (u a u†) / 2^n    = u (a / 2^n) u†
--   (p ∥ p') / 2^n    = (p / 2^n) ∥ (p' / 2^n)

seg0 :: BareAction -> Natural -> BareAction
seg0 Id                    _ = Id
seg0 (Delta k)             n = Delta (k + n)
seg0 (NegDelta k)          n = NegDelta (k + n)
seg0 (DeltaPow k j)        n = DeltaPow (k + n) j
seg0 (DeltaParam k a)      n = DeltaParam (k + n) a
seg0 (NegDeltaParam k a)   n = NegDeltaParam (k + n) a
seg0 (a :⊕: b)             n = seg0 a n :⊕: seg0 b n

seg1 :: Action -> Natural -> Action
seg1 (Bare d)    n = Bare (seg0 d n)
seg1 (Conj u a)  n = Conj u (seg1 a n)

segProc :: Proc -> Natural -> Proc
segProc (Atom a)  n = Atom (seg1 a n)
segProc (p :∥: q) n = segProc p n :∥: segProc q n

--------------------------------------------------------------------------------
-- Schedule of a process
--
--   sch(a)     = [a]
--   sch(p ∥ q) = sch(p) ++ sch(q)

sch :: Proc -> Schedule
sch (Atom a)  = [a]
sch (p :∥: q) = sch p ++ sch q

--------------------------------------------------------------------------------
-- Schedule halving
--
--   [a_1, ..., a_m] / 2^n = [a_1 / 2^n, ..., a_m / 2^n]

segSchedule :: Schedule -> Natural -> Schedule
segSchedule acts n = map (`seg1` n) acts

--------------------------------------------------------------------------------
-- Canonical schedule
--
-- sch(p) / 2^n  =  sch(p / 2^n)   (by Definitions 3 and 4 / 2 and 3)

canonicalSchedule :: Proc -> Natural -> Schedule
canonicalSchedule p n = segSchedule (sch p) n

--------------------------------------------------------------------------------
-- Schedule families
--
-- Rule (R1) of Hamilton's equational theory: any permutation of a canonical schedule is
-- an equivalent schedule. This gives the brute-force family of O(m!)
-- schedules for a process with m parallel branches.

schedulePerms :: Proc -> Natural -> [Schedule]
schedulePerms proc n = permutations (canonicalSchedule proc n)

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
