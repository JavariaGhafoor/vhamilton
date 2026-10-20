-- | Minimal VHamilton examples.
--
-- These mirror the small Hamilton implementation examples but isolate the new
-- feature: a single process term denotes a family selected by a valuation.
module ExamplesVHamilton
  ( zV
  , xV
  , zxV
  , zzV
  , instantiateAndCompare
  ) where

import Numeric.Natural (Natural)

import Syntax
  ( Valuation
  , Clifford(..)
  , Proc(..), (∥)
  , BareAction(..), (⊕), Action(..)
  )
import ExamplesChem (pauliParamProc)
import Compiler (compileTrotterWith)
import Denot (compareToIdealWith)

act :: BareAction -> Proc
act = Atom . Bare

conj :: Clifford -> BareAction -> Proc
conj u d = Atom (Conj u (Bare d))

-- | alpha Z.
zV :: String -> Proc
zV alpha =
  act (DeltaParam 0 alpha ⊕ NegDeltaParam 0 alpha)

-- | alpha X.
xV :: String -> Proc
xV alpha =
  conj H (DeltaParam 0 alpha ⊕ NegDeltaParam 0 alpha)

-- | alpha Z ∥ beta X.  This is one VHamilton process family; the parallel
-- structure is still Hamilton's additive composition.
zxV :: Proc
zxV = zV "alpha" ∥ xV "beta"

-- | gamma ZZ on two qubits.
zzV :: Proc
zzV = pauliParamProc "ZZ" 0 "gamma"

-- | Close a VHamilton family with rho, compile it, and compare to its selected
-- ideal unitary.
instantiateAndCompare :: Valuation -> Natural -> Proc -> IO Double
instantiateAndCompare rho n p = do
  let qc = compileTrotterWith rho n p
  fst <$> compareToIdealWith rho p qc
