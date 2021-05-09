import numpy as np
from qutip import qeye, sigmax, sigmay, sigmaz, tensor, create, destroy, expect

def Is(i, levels=2): return [qeye(levels) for j in range(0, i)]

def Sx(N, i): return tensor(Is(i) + [sigmax()] + Is(N - i - 1))
def Sy(N, i): return tensor(Is(i) + [sigmay()] + Is(N - i - 1))
def Sz(N, i): return tensor(Is(i) + [sigmaz()] + Is(N - i - 1))

def I(N): return Sz(N, 0)*Sz(N, 0)

def osum(lst, d=None):
    if len(lst) == 0:
        return d
    p = lst[0]
    for U in lst[1:]:
        p += U
    return p

def oprd(lst, d=None):
    if len(lst) == 0:
        return d
    p = lst[0]
    for U in lst[1:]:
        p = p*U
    return p


def gl(L, n, start=0, Opers=None):
    Sa, Sb, Sc = Sz, Sy, Sx
    if Opers is not None:
        Sa, Sb, Sc = Opers
    return oprd([Sa(L, j) for j in range(start, n)], d=I(L))*Sc(L, n)


def gr(L, n, start=0, Opers=None):
    Sa, Sb, Sc = Sz, Sy, Sx
    if Opers is not None:
        Sa, Sb, Sc = Opers
    return oprd([Sa(L, j) for j in range(start, n)], d=I(L))*Sb(L, n)


def f_(L, g1, i, g2, j, start=0, Opers=None):
    ga = g1(L, i, start=start, Opers=Opers)
    gb = g2(L, j, start=start, Opers=Opers)
    return 0.5*(ga+1j*gb)


def fd(L, g1, i, g2, j, start=0, Opers=None):
    ga = g1(L, i, start=start, Opers=Opers)
    gb = g2(L, j, start=start, Opers=Opers)
    return 0.5*(ga-1j*gb)


def a_(L, n, start=0, Opers=None):
    return f_(L, gl, n, gr, n, start=start, Opers=Opers)


def ad(L, n, start=0, Opers=None):
    return fd(L, gl, n, gr, n, start=start, Opers=Opers)


def N(L, gi, i, gj, j, start=0, Opers=None):
    _fd = fd(L, gi, i, gj, j, start=start, Opers=Opers)
    _f_ = f_(L, gi, i, gj, j, start=start, Opers=Opers)
    return _fd*_f_


def commutator(A, B):
    return A*B - B*A


def anticommutator(A, B):
    return A*B + B*A


def prepareMeasurementOperators(L, start=0, Opers=None):
    indices = []
    measure_operators = []
    for i in range(L):
        for si, gi in [('l', gl), ('r', gr)]:
            for j in range(L):
                for sj, gj in [('l', gl), ('r', gr)]:
                    label = str(i) + si + str(j) + sj
                    index = (i, si, j, sj, label)
                    operator = N(L, gi, i, gj, j, start=0, Opers=None)
                    indices.append(index)
                    measure_operators.append(operator)
    return indices, measure_operators


# function that calculates probability of measuring each pairing
def measure(psi, measure):
    ps = []
    for operator in measure:
        p = np.abs(expect(operator, psi))**2.
        ps.append(p)
    return ps


def Uij(L, g1, i, g2, j, start=0, Opers=None):
    ga = g1(L, i, start=start, Opers=Opers)
    gb = g2(L, j, start=start, Opers=Opers)
    return ((np.pi/4.)*ga*gb).expm()
