import pytest
import itertools

from qm import Sx, Sy, Sz, I, osum
from qm import f_, fd, a_, ad, gl, gr
from qm import commutator, anticommutator


def delta(ga, ia, gb, ib):
    if ga == gb and ia == ib:
        return 1
    return 0


def makeMajoranaParams(majorana_params_pre):
    majorana_params = []
    for N, maj in majorana_params_pre:
        for j in range(N):
            for i in range(j, N):
                midx = j, i
                param = N, maj, midx
                majorana_params.append(param)
    return majorana_params


def makeQuasifermionsParams(quasifermion_params_pre):
    quasifermion_params = []
    for N, maj1, maj2 in quasifermion_params_pre:
        (_, lga1), (_, lgb1) = maj1
        (_, lga2), (_, lgb2) = maj2
        for j1 in range(N):
            for j2 in range(j1, N):
                a1 = '{0}{1}'.format(j1, lga1)
                b1 = '{0}{1}'.format(j2, lgb1)
                if a1 == b1:
                    continue
                for i1 in range(N):
                    for i2 in range(i1, N):
                        # we only check fermion properties
                        # if no Majoranas are in common
                        a2 = '{0}{1}'.format(i1, lga2)
                        b2 = '{0}{1}'.format(i2, lgb2)
                        if a2 == b2:
                            continue
                        '''
                        variables = [a1, b1, a2, b2]
                        if len(set(variables)) == len(variables):
                            print(set(variables))
                            print(variables)
                            print('\n\n\n')
                        '''
                        midx1 = j1, j2
                        midx2 = i1, i2
                        params = N, maj1, midx1, maj2, midx2
                        quasifermion_params.append(params)
    return quasifermion_params


def fermionTestName(param):
    N, jw = param
    (_, la), (_, lb), (_, lc) = jw
    return 'N={0},JW={1}'.format(str(N), la + lb + lc)


def majoranaTestName(param):
    N, maj, midx = param
    la, lb, lc = 'Z','Y','X'
    (_, lga), (_, lgb) = maj
    iga, igb = midx
    return 'N={0},JW={1},(ga,gb)=({2}{3},{4}{5})'.format(
        str(N),
        la + lb + lc,
        lga,
        str(iga),
        lgb,
        str(igb))


def quasifermionTestName(param):
    N, maj1, midx1, maj2, midx2 = param
    la, lb, lc = 'Z','Y','X'
    (_, lga1), (_, lgb1) = maj1
    iga1, igb1 = midx1
    (_, lga2), (_, lgb2) = maj2
    iga2, igb2 = midx2
    return 'N={0},JW={1},f1(ga,gb)=({2}{3},{4}{5}),f2(ga,gb)=({6}{7},{8}{9})'.format(
        str(N),
        la + lb + lc,
        lga1,
        str(iga1),
        lgb1,
        str(igb1),
        lga2,
        str(iga2),
        lgb2,
        str(igb2))


Ns = [2, 3, 4]
jws = list(itertools.permutations([(Sx, 'X'), (Sy, 'Y'), (Sz, 'Z')]))
majs = list(itertools.product([(gl, 'l'), (gr, 'r')], repeat=2))

fermion_params = list(itertools.product(Ns, jws))
fermion_names = [fermionTestName(param) for param in fermion_params]

majorana_params_pre = list(itertools.product(Ns, majs))
majorana_params = makeMajoranaParams(majorana_params_pre)
majorana_names = [majoranaTestName(param) for param in majorana_params]

quasifermion_params_pre = list(itertools.product([2], majs, majs))
quasifermion_params = makeQuasifermionsParams(quasifermion_params_pre)
quasifermion_names = [quasifermionTestName(param) for param in quasifermion_params]


@pytest.mark.parametrize('N,jw', fermion_params, ids=fermion_names)
def testFermions(N, jw):
    (Sa, _), (Sb, _), (Sc, _) = jw
    Opers = Sa, Sb, Sc
    zero = 0.*I(N)
    # test all the pairs
    for n in range(N):
        a_n = a_(N, n, Opers=Opers)
        adn = ad(N, n, Opers=Opers)
        for np in range(N):
            a_np = a_(N, np, Opers=Opers)
            adnp = ad(N, np, Opers=Opers)
            assert anticommutator(a_n, a_np) == zero
            if n == np:
                assert anticommutator(a_n, adnp) == I(N)
            else:
                assert anticommutator(a_n, adnp) == zero
        assert a_n*a_n == zero
        assert adn*adn == zero


@pytest.mark.parametrize('N,maj,idxs', majorana_params, ids=majorana_names)
def testMajoranas(N,maj,idxs):
    Sa, Sb, Sc = Sz, Sy, Sx
    (g1, s1), (g2, s2) = maj
    i1, i2 = idxs
    Opers = Sa, Sb, Sc
    zero = 0.*I(N)
    # test Majorana properties
    ga = g1(N, i1, Opers=Opers)
    gb = g2(N, i2, Opers=Opers)
    if i1 == i2 and s1 == s2:
        assert anticommutator(ga, gb) == 2*I(N)
    else:
        assert anticommutator(ga, gb) == zero
    assert ga*ga == I(N)
    assert gb*gb == I(N)


@pytest.mark.parametrize('N,maj1,idxs1,maj2,idxs2', quasifermion_params, ids=quasifermion_names)
@pytest.mark.skip(reason="something is wrong with this test")
def testQuasifermions(N,maj1,idxs1,maj2,idxs2):
    (g11, s11), (g12, s12) = maj1
    i11, i12 = idxs1
    (g21, s21), (g22, s22) = maj2
    i21, i22 = idxs2
    Opers = Sz, Sy, Sx
    # build quasifermion operators
    fa_ = f_(N, g11, i11, g12, i12, Opers=Opers)
    fad = fd(N, g11, i11, g12, i12, Opers=Opers)
    fb_ = f_(N, g21, i21, g22, i22, Opers=Opers)
    fbd = fd(N, g21, i21, g22, i22, Opers=Opers)
    # test quasifermions
    # assert anticommutator(fa_, fb_) == zero
    daa = delta(g11, i11, g21, i21)
    dab = delta(g11, i11, g22, i22)
    dba = delta(g12, i12, g21, i21)
    dbb = delta(g12, i12, g22, i22)
    assert commutator(fa_, fb_) == 0.5*(daa - dbb + 1j*dab + 1j*dba)*I(N)
    assert commutator(fa_, fbd) == 0.5*(daa + dbb - 1j*dab + 1j*dba)*I(N)
