import numpy as np

from qutip import basis, tensor, Options, mesolve, qeye

from qm import Sz, Sy, Sx
from qm import gl, gr, a_, ad, osum

from qm import Uij

from plotting import plotExpectations

# glglglglglglglglgggg
url = 'https://mareknarozniak.com/2021/06/09/majorana-qubits/'

# length of the single 1D topological semiconductor
L = 2

# total system size consisting of two chains
Ltot = 2*L

# Jordan-Wigner choice to have |0> as the vacuum
Opers = Sz, Sy, Sx

# produce a vacuum state for single 1DTS of length L
vac = basis(2, 0)
vacuum = tensor([vac for j in range(L)])

# produce a completely filled state
filled = vacuum
for n in range(L):
    filled = ad(L, n, Opers=Opers)*filled

# produce the logical state of the form |0>
psi0L = filled
for i, n in enumerate(reversed(range(L-1))):
    U = Uij(L, gl, n, gl, n+1, Opers=Opers)
    psi0L = U*psi0L

# produce the logical state of the form |1>
psi1L = a_(L, L-1, Opers=Opers)*filled
for i, n in enumerate(reversed(range(L-1))):
    U = Uij(L, gl, n, gl, n+1, Opers=Opers)
    psi1L = U*psi1L

# make an initial state being a two-wire state of the logical form |01>
psi01L = tensor(psi0L, psi1L)

# prepare the simulation
tau = 2*np.pi
resolution = 60
times = np.linspace(0., tau, resolution)
# opts = Options(store_final_state=True)

# inner inter-wire braid Hamiltonian operator
H = 0.25*1j*gr(Ltot, L-1, Opers=Opers)*gl(Ltot, L, Opers=Opers)

# produce a measurement operators that measure the two nanowire system in the logical Z-basis
lSz = psi0L*psi0L.dag() - psi1L*psi1L.dag()
lSz1 = tensor([lSz] + [qeye(2) for i in range(L)])
lSz2 = tensor([qeye(2) for i in range(L)] + [lSz])

# here perform the time evolution and get the expectation values
wires = mesolve(H, psi01L, times, e_ops=[lSz1, lSz2]).expect

# as for the qubit system, we make similar time evolution just to be able to compare them
psi01 = tensor(basis(2, 0), basis(2, 1))

# qubit Hamiltonian
H = 0.25*Sx(2, 0)*Sx(2, 1)

# simulate qubits
qubits = mesolve(H, psi01, times, e_ops=[Sz(2, 0), Sz(2, 1)]).expect

# plot the results
plotExpectations(
    times,
    list(qubits) + list(wires),
    ['$<\sigma^z_1>$', '$<\sigma^z_2>$', '$<S^z_1>$', '$<S^z_2>$'],
    ['tab:blue', 'tab:orange', 'black', 'white'],
    ['-', '-', '--', '--'],
    [3., 3., 1., 1.],
    r'Comparing 1DTS $i\frac{1}{4}\gamma_{(L,r)}\gamma_{(L+1,l)}$ to qubits $\sqrt{X_1 X_2}$',
    'plots/majorana-qubits/result.png', axvlines=[np.pi])
