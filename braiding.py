import numpy as np

from qutip import basis, tensor, Options, mesolve

from qm import Sz, Sy, Sx
from qm import gl, gr, a_, ad, osum

from qm import prepareMeasurementOperators, measure

from qm import Uij

from plotting import plotMeasurement

# huehuehue
url = 'https://mareknarozniak.com/2021/05/09/braiding/'

# length of the 1D topological semiconductor
L = 3

# Jordan-Wigner choice to have |0> as the vacuum
Opers = Sz, Sy, Sx

# produce a vacuum state for single 1DTS of length L
vac = basis(2, 0)
vacuum = tensor([vac for j in range(L)])

# produce a completely filled state
filled = vacuum
for n in range(L):
    filled = ad(L, n, Opers=Opers)*filled

# prepare the measurement operators
indices, measure_operators = prepareMeasurementOperators(L, Opers=Opers)

# reproduce the effect of the topological phase transition
psi = filled
for i, n in enumerate(reversed(range(L-1))):
    U = Uij(L, gl, n, gl, n+1, Opers=Opers)
    psi = U*psi

filename = 'plots/braiding/psi0L'.format(str(i+1))
title = r'$\left\vert 0_L \right\rangle$'.format(str(i+1))
fidelities = measure(psi, measure_operators)
plotMeasurement(L, indices, fidelities, filename=filename, title=title)

# what is in |psi> is our logical zero,
psi0L = psi

# we can produce logical one state by delocalizing
# vacuum to the edges
psi = a_(L, L-1, Opers=Opers)*filled
for i, n in enumerate(reversed(range(L-1))):
    U = Uij(L, gl, n, gl, n+1, Opers=Opers)
    psi = U*psi

filename = 'plots/braiding/psi1L'.format(str(i+1))
title = r'$\left\vert 1_L \right\rangle$'.format(str(i+1))
fidelities = measure(psi, measure_operators)
plotMeasurement(L, indices, fidelities, filename=filename, title=title)

# and that is the logical one state
psi1L = psi

# once we have logical zero and one we can also define
# the eigenstates of sigma X and sigma Y operators
psip = (1/np.sqrt(2))*(psi0L + psi1L)
psim = (1/np.sqrt(2))*(psi0L - psi1L)
psipj = (1/np.sqrt(2))*(psi0L + 1j*psi1L)
psimj = (1/np.sqrt(2))*(psi0L - 1j*psi1L)

# now we have topological regime produced
# let us copy it to have two nanowires next to each other
psi0 = tensor(psi0, psi0)

# and also make a |00> qubit state for the reference
psiq0 = tensor(basis(2, 0), basis(2, 0))
