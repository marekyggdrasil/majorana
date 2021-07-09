import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D


def plotMeasurement(L, indices, measurements, s=0.25, filename=None, title=None, url=None):
    fig = plt.figure()
    axes = fig.add_subplot(111, projection='3d')
    x, y, z = [], [], []
    dx, dy, dz, c = [], [], [], []
    l = 1. - s
    cmap = plt.get_cmap('plasma')
    column_names = []
    for j in range(L):
        name = '{0}L'.format(str(j+1))
        column_names.append(name)
        name = '{0}R'.format(str(j+1))
        column_names.append(name)
    for index, measurement in zip(indices, measurements):
        i, si, j, sj, label = index
        idx_i = 2*i
        if si == 'r':
            idx_i += 1
        idx_j = 2*j
        if sj == 'r':
            idx_j += 1
        x += [idx_i]
        y += [idx_j]
        z += [0]
        dx += [l]
        dy += [l]
        dz += [measurement]
        c += [cmap(measurement)]
    ticks = np.arange(0, 2*L, 1)
    axes.set_zlim(0., 1.)
    axes.set_xticks(ticks)
    axes.set_xticklabels(column_names)
    axes.set_yticks(ticks)
    axes.set_yticklabels(column_names)
    axes.set_xlabel('source')
    axes.set_ylabel('target')
    axes.set_zlabel('$P$')
    axes.bar3d(x, y, z, dx, dy, dz, color=c, zsort='max', shade=True, edgecolor='white')
    if url is not None:
        axes.text(-0.75, 2*L+1, 0.0, url, 'y', fontsize=7)
    x_txt = axes.text(0, 0, 1.01, '$135^{\circ}$', 'x')
    y_txt = axes.text(0, 2*L-0.2, 1.015, '$45^{\circ}$', 'y')
    x_title = axes.text(2*L, 2*L, 1.2, title, 'y', fontsize=20)
    y_title = axes.text(2*L, 0, 1.25, title, 'x', fontsize=20)
    axes.view_init(elev=45, azim=45)
    if filename is not None:
        x_txt.set_visible(False)
        y_txt.set_visible(True)
        x_title.set_visible(False)
        y_title.set_visible(True)
        fig.savefig(filename + '_45.png', transparent=True, bbox_inches='tight', pad_inches=0)
    axes.view_init(elev=45, azim=135)
    if filename is not None:
        x_txt.set_visible(True)
        y_txt.set_visible(False)
        x_title.set_visible(True)
        y_title.set_visible(False)
        fig.savefig(filename + '_135.png', transparent=True, bbox_inches='tight', pad_inches=0)


def plotExpectations(times, expects, labels, colors, styles, linewidths, title, filename, axvlines=[]):
    fig, ax = plt.subplots(1, 1, constrained_layout=True)
    ax.set_title(title)
    ax.set_xlabel('t')
    ax.grid()
    for x in axvlines:
        ax.axvline(x=x, color='black', linestyle='--')
    for measurement, label, color, style, width in zip(expects, labels, colors, styles, linewidths):
        ax.plot(times, measurement, label=label, color=color, linestyle=style, linewidth=width)
    ax.legend()
    fig.savefig(filename, transparent=True)


# from: https://matplotlib.org/stable/gallery/lines_bars_and_markers/barchart.html#sphx-glr-gallery-lines-bars-and-markers-barchart-py
def plotTeleportationOutcomes(outcomes, corrs, labels, title, url=None, filename=None):
    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots(1, 1, constrained_layout=True)
    rects1 = ax.bar(x - width/2, outcomes, width, label='not corrected')
    rects2 = ax.bar(x + width/2, corrs, width, label='corrected')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('fidelity')
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()

    ax.bar_label(rects1, padding=3)
    ax.bar_label(rects2, padding=3)

    if url is not None:
        ax.text(0, -0.1, url, fontsize=7)

    fig.tight_layout()
    if filename is not None:
        fig.savefig(filename, transparent=True)
