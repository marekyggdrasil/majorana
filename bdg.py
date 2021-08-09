import numpy as np

from qm import Hbdg
from plotting import plotHbdgSpectrum


L = 20
mu_max = 4.0
t = 1.0
delta = 1.0

res = 20

mus = np.linspace(0., mu_max, res)
spectrum = np.zeros((res, 2*L))

for j, mu in enumerate(mus):
    H = Hbdg(L, mu, t, delta)
    eigenvalues, _ = np.linalg.eigh(H)
    spectrum[j, :] = eigenvalues

url = 'https://mareknarozniak.com/2021/08/09/bdg/'
title = 'Topological gap in Bogoliubov-de-Gennes formalism'
plotHbdgSpectrum(L, mus, spectrum, title, mark=2.0, url=url, url_x=0.0, url_y=4.5, filename='plots/bdg/result.png')
