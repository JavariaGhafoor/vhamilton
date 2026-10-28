------------------------------------------------------------------------
-- VHamilton: open clocks over Hamilton, with a symbolic Hamiltonian
-- semantics in the Pauli basis, and a machine-checked equivalence of
-- two variational programs for every valuation and every clock subscript k.
--
-- Checked with Agda 2.6.3 and agda-stdlib 1.7.3:
--   agda -i . -i /usr/share/agda-stdlib VHamilton.agda
-- Written to avoid the stdlib 1.x / 2.x difference in Data.Vec.replicate.
------------------------------------------------------------------------

module VHamilton where

open import Data.Nat     using (ℕ; zero; suc; _+_; _^_)
open import Data.Fin     using (Fin; zero; suc)
open import Data.Vec     using (Vec; []; _∷_; take; drop; _[_]≔_; toList; tabulate)
                         renaming (_++_ to _++ᵥ_; map to mapᵥ; zipWith to zipWithᵥ)
open import Data.List    using (List; []; _∷_; _++_; map; foldr; concat; zipWith)
open import Data.Product using (_×_; _,_)
open import Data.Bool    using (Bool; true; false; not; if_then_else_; _∧_)
open import Data.Integer using (+_)
open import Data.Rational using (ℚ; 0ℚ; _/_)
                         renaming (_+_ to _+ℚ_; _*_ to _*ℚ_; -_ to -ℚ_)
open import Relation.Binary.PropositionalEquality using (_≡_; refl; cong)

------------------------------------------------------------------------
-- Pauli strings, the basis in which we read Hamiltonians.

data Pauli : Set where
  𝟙 X Y Z : Pauli

PString : ℕ → Set
PString = Vec Pauli

eqP : Pauli → Pauli → Bool
eqP 𝟙 𝟙 = true
eqP X X = true
eqP Y Y = true
eqP Z Z = true
eqP _ _ = false

eqS : ∀ {n} → PString n → PString n → Bool
eqS []       []       = true
eqS (a ∷ as) (b ∷ bs) = eqP a b ∧ eqS as bs

------------------------------------------------------------------------
-- Transformations (Hamilton's Clifford layer), acting on Pauli strings
-- by conjugation P ↦ u P u†, which yields ±P'.  The sign is a Bool
-- (true = +).  Tables were checked against the matrices.

data Cliff : ℕ → Set where
  I₁ H S   : Cliff 1
  CX SWAP  : Cliff 2
  _·_      : ∀ {n}   → Cliff n → Cliff n → Cliff n
  _⊗_      : ∀ {m n} → Cliff m → Cliff n → Cliff (m + n)

sgn : Bool → Bool → Bool
sgn s t = if s then t else not t

cjH : Pauli → Bool × Pauli
cjH 𝟙 = true , 𝟙
cjH X = true , Z
cjH Y = false , Y
cjH Z = true , X

cjS : Pauli → Bool × Pauli
cjS 𝟙 = true , 𝟙
cjS X = true , Y
cjS Y = false , X
cjS Z = true , Z

-- CX with the first qubit as control.
cjCX : Pauli → Pauli → Bool × PString 2
cjCX 𝟙 𝟙 = true  , 𝟙 ∷ 𝟙 ∷ []
cjCX 𝟙 X = true  , 𝟙 ∷ X ∷ []
cjCX 𝟙 Y = true  , Z ∷ Y ∷ []
cjCX 𝟙 Z = true  , Z ∷ Z ∷ []
cjCX X 𝟙 = true  , X ∷ X ∷ []
cjCX X X = true  , X ∷ 𝟙 ∷ []
cjCX X Y = true  , Y ∷ Z ∷ []
cjCX X Z = false , Y ∷ Y ∷ []
cjCX Y 𝟙 = true  , Y ∷ X ∷ []
cjCX Y X = true  , Y ∷ 𝟙 ∷ []
cjCX Y Y = false , X ∷ Z ∷ []
cjCX Y Z = true  , X ∷ Y ∷ []
cjCX Z 𝟙 = true  , Z ∷ 𝟙 ∷ []
cjCX Z X = true  , Z ∷ X ∷ []
cjCX Z Y = true  , 𝟙 ∷ Y ∷ []
cjCX Z Z = true  , 𝟙 ∷ Z ∷ []

one : (Pauli → Bool × Pauli) → PString 1 → Bool × PString 1
one f (a ∷ []) with f a
... | s , b = s , b ∷ []

cj : ∀ {n} → Cliff n → PString n → Bool × PString n
cj I₁   P               = true , P
cj H    P               = one cjH P
cj S    P               = one cjS P
cj CX   (a ∷ b ∷ [])    = cjCX a b
cj SWAP (a ∷ b ∷ [])    = true , b ∷ a ∷ []
cj (u · v) P with cj v P
... | s , P′ with cj u P′
...   | t , P″ = sgn s t , P″
cj (_⊗_ {m} u v) P with cj u (take m P) | cj v (drop m P)
... | s , P₁ | t , P₂ = sgn s t , P₁ ++ᵥ P₂

------------------------------------------------------------------------
-- Coefficients are affine forms in the clock values: a constant plus
-- one coefficient per variable.  Variable α stands for ν(α)/2^(k α),
-- so one symbolic computation covers every valuation ν and every k.

LF : ℕ → Set
LF m = ℚ × Vec ℚ m

-- The zero vector, via tabulate: `replicate` takes its length explicitly
-- in agda-stdlib 2.x but implicitly in 1.x, so we avoid it.
zeros : ∀ {m} → Vec ℚ m
zeros = tabulate (λ _ → 0ℚ)

lf0 : ∀ {m} → LF m
lf0 = 0ℚ , zeros

_⊞_ : ∀ {m} → LF m → LF m → LF m
(c , v) ⊞ (d , w) = c +ℚ d , zipWithᵥ _+ℚ_ v w

scale : ∀ {m} → ℚ → LF m → LF m
scale q (c , v) = q *ℚ c , mapᵥ (q *ℚ_) v

neg half : ∀ {m} → LF m → LF m
neg  = scale (-ℚ (+ 1 / 1))
half = scale (+ 1 / 2)

------------------------------------------------------------------------
-- Syntax.  Scalars are fixed clocks ±δ_k or open clocks ±δ^α; each
-- variable has one k, given separately.  Bare actions are balanced
-- direct sums on qubits.

data Scalar (m : ℕ) : Set where
  zer      : Scalar m
  δ -δ     : ℕ → Scalar m
  δ^ -δ^   : Fin m → Scalar m

data Diag (m : ℕ) : ℕ → Set where
  leaf : Scalar m → Diag m 0
  _⊕_  : ∀ {n} → Diag m n → Diag m n → Diag m (suc n)

data Act (m n : ℕ) : Set where
  bare : Diag m n → Act m n
  conj : Cliff n → Act m n → Act m n      -- u a u†

data Proc (m n : ℕ) : Set where
  act : Act m n → Proc m n
  _∥_ : Proc m n → Proc m n → Proc m n

infixr 5 _⊕_
infixr 4 _∥_

------------------------------------------------------------------------
-- Hamiltonian semantics, as a list of Pauli terms.

pow½ : ℕ → ℚ
pow½ zero    = + 1 / 1
pow½ (suc k) = (+ 1 / 2) *ℚ pow½ k

unit : ∀ {m} → Fin m → ℚ → LF m
unit α q = 0ℚ , (zeros [ α ]≔ q)

valS : ∀ {m} → Scalar m → LF m
valS zer     = lf0
valS (δ k)   = pow½ k , zeros
valS (-δ k)  = -ℚ pow½ k , zeros
valS (δ^ α)  = unit α (+ 1 / 1)
valS (-δ^ α) = unit α (-ℚ (+ 1 / 1))

Ham : ℕ → ℕ → Set
Ham m n = List (LF m × PString n)

-- d₀ ⊕ d₁ = |0⟩⟨0| ⊗ d₀ + |1⟩⟨1| ⊗ d₁ = 𝟙 ⊗ (d₀+d₁)/2 + Z ⊗ (d₀−d₁)/2
split : ∀ {m n} → LF m × PString n → LF m × PString n → Ham m (suc n)
split (c₀ , P) (c₁ , _) = (half (c₀ ⊞ c₁) , 𝟙 ∷ P) ∷ (half (c₀ ⊞ neg c₁) , Z ∷ P) ∷ []

diag : ∀ {m n} → Diag m n → Ham m n
diag (leaf s)  = (valS s , []) ∷ []
diag (d₀ ⊕ d₁) = concat (zipWith split (diag d₀) (diag d₁))

conjT : ∀ {m n} → Cliff n → LF m × PString n → LF m × PString n
conjT u (c , P) with cj u P
... | true  , P′ = c , P′
... | false , P′ = neg c , P′

⟦_⟧A : ∀ {m n} → Act m n → Ham m n
⟦ bare d   ⟧A = diag d
⟦ conj u a ⟧A = map (conjT u) ⟦ a ⟧A

⟦_⟧ : ∀ {m n} → Proc m n → Ham m n
⟦ act a ⟧ = ⟦ a ⟧A
⟦ p ∥ q ⟧ = ⟦ p ⟧ ++ ⟦ q ⟧

-- Coefficient of one Pauli string.
coef : ∀ {m n} → Ham m n → PString n → LF m
coef h Q = foldr step lf0 h
  where step : _ → LF _ → LF _
        step (c , P) acc = if eqS P Q then c ⊞ acc else acc

------------------------------------------------------------------------
-- Valuations.  A valuation fills each open clock with a natural number;
-- k α is the subscript of the clock δ_k^α (one per variable).

Val : ℕ → Set
Val m = Fin m → ℕ

evalLF : ∀ {m} → (Fin m → ℕ) → Val m → LF m → ℚ
evalLF k ν (c , v) =
  foldr _+ℚ_ c (toList (zipWithᵥ _*ℚ_ v (tabulate λ α → (+ ν α / 1) *ℚ pow½ (k α))))

-- ⟦ p ⟧^ν : the coefficient of each Pauli string under a valuation.
⟦_⟧^ : ∀ {m n} → Proc m n → (Fin m → ℕ) → Val m → PString n → ℚ
⟦ p ⟧^ k ν Q = evalLF k ν (coef ⟦ p ⟧ Q)

-- Hamiltonian equivalence: equal coefficients, as affine forms.
_≃_ : ∀ {m n} → Proc m n → Proc m n → Set
p ≃ q = ∀ Q → coef ⟦ p ⟧ Q ≡ coef ⟦ q ⟧ Q

-- It implies equality under every valuation ν and every k.
sound : ∀ {m n} {p q : Proc m n} → p ≃ q →
        ∀ k ν Q → ⟦ p ⟧^ k ν Q ≡ ⟦ q ⟧^ k ν Q
sound eq k ν Q = cong (evalLF k ν) (eq Q)

------------------------------------------------------------------------
-- Example: one clock γ on two qubits.

γ : Fin 1
γ = zero

-- ZZ(γ) = δ^γ ⊕ −δ^γ ⊕ −δ^γ ⊕ δ^γ, i.e. (ν(γ)/2^k) Z⊗Z.
ZZ : Proc 1 2
ZZ = act (bare ((leaf (δ^ γ) ⊕ leaf (-δ^ γ)) ⊕ (leaf (-δ^ γ) ⊕ leaf (δ^ γ))))

-- I⊗Z(γ) = δ^γ ⊕ −δ^γ ⊕ δ^γ ⊕ −δ^γ, conjugated by CX.
IZ : Act 1 2
IZ = bare ((leaf (δ^ γ) ⊕ leaf (-δ^ γ)) ⊕ (leaf (δ^ γ) ⊕ leaf (-δ^ γ)))

CX-IZ-CX : Proc 1 2
CX-IZ-CX = act (conj CX IZ)

-- The two programs are equivalent: every Pauli coefficient agrees.
-- Split on the first Pauli only; Agda computes both coefficients of the
-- four strings in each clause, and refl closes all 16.
ZZ≃CX-IZ-CX : ZZ ≃ CX-IZ-CX
ZZ≃CX-IZ-CX (𝟙 ∷ _ ∷ []) = refl
ZZ≃CX-IZ-CX (X ∷ _ ∷ []) = refl
ZZ≃CX-IZ-CX (Y ∷ _ ∷ []) = refl
ZZ≃CX-IZ-CX (Z ∷ _ ∷ []) = refl

-- Hence they denote the same Hamiltonian for every valuation ν and every k.
ZZ-equiv : ∀ k ν Q → ⟦ ZZ ⟧^ k ν Q ≡ ⟦ CX-IZ-CX ⟧^ k ν Q
ZZ-equiv = sound {p = ZZ} {q = CX-IZ-CX} ZZ≃CX-IZ-CX

-- Sanity check: the coefficient of Z⊗Z is exactly the symbol γ (1·γ).
ZZ-coef : coef ⟦ ZZ ⟧ (Z ∷ Z ∷ []) ≡ unit γ (+ 1 / 1)
ZZ-coef = refl
