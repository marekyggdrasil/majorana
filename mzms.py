import numpy as np

from qutip import basis, tensor, Options, mesolve, expect

from qm import Sz, Sy, Sx
from qm import gl, gr, osum
from qm import a_ as op_a_
from qm import ad as op_ad
from qm import f_ as op_f_
from qm import fd as op_fd

from plotting import plotMeasurement

# huehuehue
url = 'https://mareknarozniak.com/2021/04/09/mzms/'

# length of the 1D topological semiconductor
L = 3

# Jordan-Wigner choice to have |0> as the vacuum
Opers = Sz, Sy, Sx

# prepare operators for more readable notation
a_ = lambda n: op_a_(L, n, Opers=Opers)
ad = lambda n: op_ad(L, n, Opers=Opers)
f_ = lambda gi, i, gj, j: op_f_(L, gi, i, gj, j, Opers=Opers)
fd = lambda gi, i, gj, j: op_fd(L, gi, i, gj, j, Opers=Opers)
N = lambda gi, i, gj, j: fd(gi, i, gj, j)*f_(gi, i, gj, j)

# produce a vacuum state

vac = basis(2, 0)
vacuum = tensor([vac for j in range(L)])

# produce a completely filled state
filled = vacuum
for n in range(L):
    filled = ad(n)*filled

# produce a Kitaev chain Hamiltonian
H_trivial = osum([ad(n)*a_(n) for n in range(L)])
Hhop = lambda n: ad(n+1)*a_(n)+ad(n)*a_(n+1)
Hbcs = lambda n: a_(n)*a_(n+1)+ad(n+1)*ad(n)
H_hopping = osum([ad(n+1)*a_(n)+ad(n)*a_(n+1) for n in range(L-1)])
H_bcs = osum([a_(n)*a_(n+1)+ad(n+1)*ad(n) for n in range(L-1)])
H = lambda mu, t, d: -0.5*mu*H_trivial-t*H_hopping+d*H_bcs

# check the energy of the trivial regime
energy_trivial = expect(H(2., 0., 0.), filled)
print('Energy of the trivial regime is', energy_trivial)

# prepare the measurement operators
indices = []
measure_operators = []
for i in range(L):
    for si, gi in [('l', gl), ('r', gr)]:
        for j in range(L):
            for sj, gj in [('l', gl), ('r', gr)]:
                label = str(i) + si + str(j) + sj
                index = (i, si, j, sj, label)
                operator = N(gi, i, gj, j)
                indices.append(index)
                measure_operators.append(operator)

# function that calculates probability of measuring each pairing
def measure(psi, measure):
    ps = []
    for operator in measure:
        p = np.abs(expect(operator, psi))**2.
        ps.append(p)
    return ps

# prepare the simulation
tau = 10.*np.pi
resolution = 300
times = np.linspace(0., tau, resolution)
opts = Options(store_final_state=True)

psi0 = filled
fidelities = measure(psi0, measure_operators)
plotMeasurement(L, indices, fidelities, filename='plots/mzms/psi0', title=r'$\left\vert \psi_0 \right\rangle$')

H0 = -0.5*(2*ad(0)*a_(0)+2*ad(1)*a_(1)+2*ad(2)*a_(2)-3)
Hf = -0.5*(2*ad(0)*a_(0)+2*ad(1)*a_(1)-1)
Ht = lambda t, args: (1-t/tau)*H0 + (t/tau)*Hf
psi1 = mesolve(Ht, psi0, times, options=opts).final_state
fidelities = measure(psi1, measure_operators)
plotMeasurement(L, indices, fidelities, filename='plots/mzms/psi1', title=r'$\left\vert \psi_1 \right\rangle$')

H0 = -0.5*(2*ad(0)*a_(0)+2*ad(1)*a_(1)-2)
Hf = -0.5*(2*ad(0)*a_(0)-1)+Hhop(1)-Hbcs(1)
Ht = lambda t, args: (1-t/tau)*H0 + (t/tau)*Hf
psi2 = mesolve(Ht, psi1, times, options=opts).final_state
fidelities = measure(psi2, measure_operators)
plotMeasurement(L, indices, fidelities, filename='plots/mzms/psi2', title=r'$\left\vert \psi_2 \right\rangle$', url=url)

H0 = -0.5*(2*ad(0)*a_(0)-1)+Hhop(1)-Hbcs(1)
Hf = Hhop(0)-Hbcs(0)+Hhop(1)-Hbcs(1)
Ht = lambda t, args: (1-t/tau)*H0 + (t/tau)*Hf
psi3 = mesolve(Ht, psi2, times, options=opts).final_state
fidelities = measure(psi3, measure_operators)
plotMeasurement(L, indices, fidelities, filename='plots/mzms/psi3', title=r'$\left\vert \psi_3 \right\rangle$')
