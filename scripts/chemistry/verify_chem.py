"""
Two DISTINCT error sources, measured separately and honestly:
 (A) dyadic representation error: || exp(-2pi i H_true) - exp(-2pi i H_dyadic(k)) ||
     -- depends only on k, not n. The cost of rounding real coeffs to level k.
 (B) Trotter error at level n, against the dyadic Hamiltonian it actually targets.
The TOTAL error against the true dynamics is bounded by (A)+(B); we also report it
directly.
"""
import json, subprocess
import numpy as np
from math import pi
from itertools import product
from scipy.linalg import expm

I=np.eye(2);X=np.array([[0,1],[1,0]],dtype=complex)
Y=np.array([[0,-1j],[1j,0]]);Z=np.array([[1,0],[0,-1]],dtype=complex)
P={'I':I,'X':X,'Y':Y,'Z':Z}
Hd=(1/np.sqrt(2))*np.array([[1,1],[1,-1]],dtype=complex);Sg=np.array([[1,0],[0,1j]],dtype=complex)
def pmat(l):
    M=np.array([[1]],dtype=complex)
    for c in l: M=np.kron(M,P[c])
    return M
def signed_exps(lbl,num):
    zmask=[c in 'XYZ' for c in lbl]; o=[]
    for (i,j) in [(0,0),(0,1),(1,0),(1,1)]:
        o.append(num if sum(b for b,m in zip([i,j],zmask) if m)%2==0 else -num)
    return o
def core_ham(lbl,num,k,n):
    return np.diag([j/2**(k+n) for j in signed_exps(lbl,num)]).astype(complex)
def basis_cliff(lbl):
    if all(c not in 'XY' for c in lbl): return None
    facs=[Hd if c=='X' else (Sg@Hd if c=='Y' else I) for c in lbl]
    U=np.array([[1]],dtype=complex)
    for f in facs: U=np.kron(U,f)
    return U
def term_U(lbl,num,k,n):
    core=core_ham(lbl,num,k,n);U=basis_cliff(lbl)
    return expm(-2j*pi*(core if U is None else U@core@U.conj().T))
def dyadic_ham(terms):
    H=np.zeros((4,4),dtype=complex)
    for t in terms:
        core=core_ham(t['pauli'],t['num'],t['k'],0);U=basis_cliff(t['pauli'])
        H += core if U is None else U@core@U.conj().T
    return H
def frob_al(Uc,Ui):
    inner=np.trace(Ui.conj().T@Uc)
    Uc=Uc/(inner/abs(inner)) if abs(inner)>1e-14 else Uc
    return float(np.linalg.norm(Uc-Ui))
def run(level):
    out=subprocess.run(["python3", "-m", "scripts.chemistry.chem_decompose", "--block", "upper", "--level", str(level)],
                       capture_output=True,text=True).stdout
    return json.loads(out)

# true scaled Hamiltonian (real coeffs), same scale as decompose uses
upper=np.array([
 [37.68,-7.537,0.7182,-0.0024],[-7.5367,7.182,-7.752,0.8407],
 [0.7182,-7.752,0.2725,-6.909],[-0.0024,0.8407,-6.909,-7.460]])
upper=(upper+upper.T)/2
d8=run(8); s=d8['scale']
H_true=s*upper
# remove identity (II) contribution from H_true too, since decompose drops it (global phase)
# (II coeff just shifts global phase; align_phase removes it, so leave as is)
U_true=expm(-2j*pi*H_true)

print("=== (A) Dyadic representation error vs level k (independent of n) ===")
for k in (3,4,5,6,8,10,12):
    terms=run(k)['terms']
    U_dy=expm(-2j*pi*dyadic_ham(terms))
    print(f"  k={k:2d}:  || U_true - U_dyadic ||_aligned = {frob_al(U_dy,U_true):.5f}")
print()
print("=== (B) Trotter error vs n, against the TRUE dynamics, at k=12 (fine coeffs) ===")
terms=run(12)['terms']
print(f"{'n':>3} {'reps':>6} {'Frob vs true':>14}")
for n in range(0,11):
    step=np.eye(4,dtype=complex)
    for t in terms: step=term_U(t['pauli'],t['num'],t['k'],n)@step
    Uc=np.linalg.matrix_power(step,2**n)
    print(f"{n:>3} {2**n:>6} {frob_al(Uc,U_true):>14.5f}")
