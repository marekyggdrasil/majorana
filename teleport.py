import numpy as np

from qutip import basis, tensor, rand_ket, qeye, expect

from qm import Sz, Sy, Sx
from qm import gl, gr, a_, ad, osum

from qm import Uij

from plotting import plotTeleportationOutcomes

# huehuehue
url = 'https://mareknarozniak.com/2021/07/09/topological-teleportation/'

# length of the 1D topological semiconductor
L = 2

# for projections, we will need an operator the ignores the state
ign = tensor([qeye(2) for j in range(L)])

# Jordan-Wigner choice to have |0> as the vacuum
Opers = Sz, Sy, Sx

# produce a vacuum state for single 1DTS of length L
vac = basis(2, 0)
vacuum = tensor([vac for j in range(L)])

# produce a completely filled state
filled = vacuum
for n in range(L):
    filled = ad(L, n, Opers=Opers)*filled

# reproduce the effect of the topological phase transition
psi = filled
for i, n in enumerate(reversed(range(L-1))):
    U = Uij(L, gl, n, gl, n+1, Opers=Opers)
    psi = U*psi

# what is in |psi> is our logical zero,
psi0L = psi

# we can produce logical one state by delocalizing
# vacuum to the edges
psi = a_(L, L-1, Opers=Opers)*filled
for i, n in enumerate(reversed(range(L-1))):
    U = Uij(L, gl, n, gl, n+1, Opers=Opers)
    psi = U*psi

# and that is the logical one state
psi1L = psi

# prepare a random state to be teleported
# as done in https://mareknarozniak.com/2020/03/22/simulating-quantum-teleportation/
rand_psi = rand_ket(2)

# and extract its amplitudes
rand_a = rand_psi.overlap(basis(2, 0))
rand_b = rand_psi.overlap(basis(2, 1))

# from that we can build a topological random state to be teleported
rand_psi_topo = rand_a*psi0L + rand_b*psi1L

# now prepare the Hilbert space consisting of five topological qubits
psi_pre = tensor(rand_psi_topo, psi0L, psi0L, psi0L, psi0L)

# state preparation, make the ancilla topological qubits in the logical |+> state
U = Uij(5*L, gr, 7, gl, 8, Opers=Opers)
psi0 = U*psi_pre

# perform the actual teleportation of the random state by applying the braids
U = Uij(5*L, gr, 1, gl, 2, Opers=Opers)
psi1 = U*psi0

U = Uij(5*L, gr, 3, gl, 4, Opers=Opers)
psi2 = U*psi1

# finally, create the projections of the outcomes
P00 = tensor(psi0L.proj(), psi0L.proj(), ign, ign, ign)
P01 = tensor(psi0L.proj(), psi1L.proj(), ign, ign, ign)
P10 = tensor(psi1L.proj(), psi0L.proj(), ign, ign, ign)
P11 = tensor(psi1L.proj(), psi1L.proj(), ign, ign, ign)

# we can use the above projections to simulate the measurement outcomes
psi_out_00 = P00*psi2
psi_out_01 = P01*psi2
psi_out_10 = P10*psi2
psi_out_11 = P11*psi2

# we prepare classical correction using the double braids
U = Uij(5*L, gl, 4, gr, 5, Opers=Opers)
Z = U*U

U = Uij(5*L, gr, 5, gl, 6, Opers=Opers)
X = U*U

# we apply the corrections according to out truth table
psi_out_00_corr = Z*psi_out_00
psi_out_01_corr = X*psi_out_01
psi_out_10_corr = Z*X*psi_out_10
psi_out_11_corr = psi_out_11

# let us trace out the third topological qubit which is the teleported state
# we do it for both corrected and not corrected cases
rho_00 = psi_out_00.ptrace([4, 5])
rho_01 = psi_out_01.ptrace([4, 5])
rho_10 = psi_out_10.ptrace([4, 5])
rho_11 = psi_out_11.ptrace([4, 5])

rho_00_corr = psi_out_00_corr.ptrace([4, 5])
rho_01_corr = psi_out_01_corr.ptrace([4, 5])
rho_10_corr = psi_out_10_corr.ptrace([4, 5])
rho_11_corr = psi_out_11_corr.ptrace([4, 5])

# normalize everything after tracing out
rho_00 = rho_00 / rho_00.norm()
rho_01 = rho_01 / rho_01.norm()
rho_10 = rho_10 / rho_10.norm()
rho_11 = rho_11 / rho_11.norm()

rho_00_corr = rho_00_corr / rho_00_corr.norm()
rho_01_corr = rho_01_corr / rho_01_corr.norm()
rho_10_corr = rho_10_corr / rho_10_corr.norm()
rho_11_corr = rho_11_corr / rho_11_corr.norm()

# check fidelities of both corrected and not corrected outcomes
fid00 = expect(rho_00, rand_psi_topo)
fid01 = expect(rho_01, rand_psi_topo)
fid10 = expect(rho_10, rand_psi_topo)
fid11 = expect(rho_11, rand_psi_topo)

fid00_corr = expect(rho_00_corr, rand_psi_topo)
fid01_corr = expect(rho_01_corr, rand_psi_topo)
fid10_corr = expect(rho_10_corr, rand_psi_topo)
fid11_corr = expect(rho_11_corr, rand_psi_topo)

# plot the results
outcomes = [fid00, fid01, fid10, fid11]
corrected_outcomes = [fid00_corr, fid01_corr, fid10_corr, fid11_corr]

outcomes = [np.round(out, decimals=2) for out in outcomes]
corrected_outcomes = [np.round(out, decimals=2) for out in corrected_outcomes]

labels = [r'$\left \vert 00_L \right \rangle$',
          r'$\left \vert 01_L \right \rangle$',
          r'$\left \vert 10_L \right \rangle$',
          r'$\left \vert 11_L \right \rangle$']

title = 'Effect of classical correction on fidelities per outcome'
filename = 'plots/topological-teleportation/outcomes.png'
plotTeleportationOutcomes(outcomes, corrected_outcomes, labels, title, url=url, filename=filename)
