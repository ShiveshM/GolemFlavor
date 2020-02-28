#! /usr/bin/env python
"""
Make triangle plot
"""

from __future__ import absolute_import, division

import sys
sys.path.extend(['.', '..'])

import numpy as np
from scipy import interpolate

import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.legend_handler import HandlerPatch
plt.style.use('./paper.mplstyle')

import ternary

from utils import fr as fr_utils
from utils import misc as misc_utils
from utils.enums import MixingScenario


OUTFOLDER = './plots/'
OUTFORMAT = ['png', 'pdf']
FONTSIZE = 23
PURPLE = [0.5019607843137255, 0.0, 0.5019607843137255]


class HandlerCircle(HandlerPatch):
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height,
                       fontsize, trans):
        r = 10
        x = r + width//2 + 10
        y = height//2 - 3

        # create 
        p = Circle(xy=(x, y), radius=r)

        # update with data from oryginal object
        self.update_prop(p, orig_handle, legend)

        # move xy to legend
        p.set_transform(trans)

        return [p]


def swap(x):
    y = x.T
    return np.vstack([y[1], y[-1], y[0]]).T


def interp(i):
    x, y, z = i.T
    p = np.linspace(0, 1, len(x))
    q = np.linspace(0, 1, 100)
    s = 0.2
    x_ = interpolate.splev(q, interpolate.splrep(p, x, s=s))
    y_ = interpolate.splev(q, interpolate.splrep(p, y, s=s))
    z_ = interpolate.splev(q, interpolate.splrep(p, z, s=s))
    return np.vstack([x_, y_, z_]).T


def setupfigure(tax, show_120=True, show_001=False, legend=True):
    # Draw Boundary and Gridlines
    tax.boundary(linewidth=2.0)
    tax.gridlines(color='grey', multiple=scale/5., linewidth=1.0, alpha=0.4,
                  ls='--')
    tax.gridlines(color='grey', multiple=scale/10., linewidth=0.5, alpha=0.4,
                  ls=':')

    # Set Axis labels and Title
    tax.left_axis_label(r"$f_{\tau}^{\oplus}$", fontsize=FONTSIZE+8,
                        offset=0.2, rotation=0)
    tax.right_axis_label(r"$f_{\mu}^{\oplus}$", fontsize=FONTSIZE+8,
                         offset=0.2, rotation=0)
    tax.bottom_axis_label(r"$f_{e}^{\oplus}$", fontsize=FONTSIZE+8,
                          position=(0.55, -0.20/2, 0.5), rotation=0)

    # Remove default Matplotlib axis
    tax.get_axes().axis('off')
    tax.clear_matplotlib_ticks()

    # Legend
    if legend:
        l_size = FONTSIZE
        legend_elements = []
        if show_120:
            legend_elements.append(
                Circle((0., 0.), 0.1, facecolor='red', alpha=0.7, edgecolor='k',
                       linewidth=2., label=r'$\left (1:2:0\right )$')
            )
        legend_elements.append(
            Circle((0., 0.), 0.1, facecolor='green', alpha=0.7, edgecolor='k',
                   linewidth=2., label=r'$\left (0:1:0\right )$')
        )
        legend_elements.append(
            Circle((0., 0.), 0.1, facecolor='blue', alpha=0.7, edgecolor='k',
                   linewidth=2., label=r'$\left (1:0:0\right )$')
        )
        if show_001:
            legend_elements.append(
                Circle((0., 0.), 0.1, facecolor=PURPLE+[0.7], edgecolor='k',
                       linewidth=2., label=r'$\left (0:0:1\right )$')
            )
        legend = plt.legend(handles=legend_elements, loc=(0.65, 0.8),
                            title='Source composition',
                            fontsize=l_size,
                            handler_map={Circle: HandlerCircle()})
        plt.setp(legend.get_title(), fontsize=l_size)
        legend.get_frame().set_linestyle('-')

    # Set ticks
    tax.ticks(axis='blr', multiple=scale/5., linewidth=1, offset=0.03,
              fontsize=FONTSIZE, tick_formats='%.1f')

    tax._redraw_labels()

    return


# Load contours
bc = np.genfromtxt('./bayes_contours.csv', delimiter=',', skip_header=2)
print bc

contour_68_upper = interp(swap(bc[:,:3]))
contour_68_lower = interp(swap(bc[:,3:6]))
contour_90_upper = interp(swap(bc[:,6:9]))
contour_90_lower = interp(swap(bc[:,9:]))

# Load data
sm_010 = np.load('./mcmc_DIM3_sfr_0_1_0_mfr_3.29685E+15_7.88411E+15_6.83344E+15_fix_scale_0.1_sigma_100000_proc.npy')
sm_100 = np.load('./mcmc_DIM3_sfr_1_0_0_mfr_1.98171E+16_6.59371E+15_9.61795E+15_fix_scale_0.1_sigma_100000_proc.npy')
bsm_010 = np.load('./mcmc_chains_DIM3_sfr_0_1_0_mfr_0_1_0_MixingScenario.T13_sigma_010.npy')
bsm_100 = np.load('./mcmc_chains_DIM3_sfr_1_0_0_mfr_1_0_0_MixingScenario.T23_sigma_010.npy')

# Figure
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(111)
ax.set_aspect('equal')

# Boundary and Gridlines
scale = 1
fig, tax = ternary.figure(ax=ax, scale=scale)
setupfigure(tax, show_120=False)

# Plot
tax.scatter(sm_010, marker='.', s=0.2, alpha=0.2, color='green')
tax.scatter(sm_100, marker='.', s=0.2, alpha=0.2, color='blue')
tax.scatter(bsm_010, marker='.', s=0.2, alpha=0.2, color='green')
tax.scatter(bsm_100, marker='.', s=0.2, alpha=0.2, color='blue')

# Contour
tax.plot(contour_68_upper, linewidth=2.3, color='grey', zorder=0, alpha=0.6)
tax.plot(contour_68_lower, linewidth=2.3, color='grey', zorder=0, alpha=0.6)
tax.plot(contour_90_upper, linewidth=2.3, color='darkgrey', zorder=0, alpha=0.6)
tax.plot(contour_90_lower, linewidth=2.3, color='darkgrey', zorder=0, alpha=0.6)

# Lines
marker_style = dict(
    linestyle=' ', color='darkorange', linewidth=1.2,
    markersize=14, zorder=0
)

p1 = (0.18301213, 0.43765598, 0.37933189)
p2 = (0, 1, 0)
divisions = 46
x_d = np.linspace(p1[0], p2[0], divisions+1)
y_d = np.linspace(p1[1], p2[1], divisions+1)
z_d = np.linspace(p1[2], p2[2], divisions+1)
for n in xrange(divisions-2):
    p = (x_d[n], y_d[n], z_d[n])
    q = (x_d[n+1], y_d[n+1], z_d[n+1])
    tax.line(p, q, marker=(3, 2, 180+46.5), **marker_style)

p1 = (0.55003613, 0.18301213, 0.26695174)
p2 = (1, 0, 0)
divisions = 35
x_d = np.linspace(p1[0], p2[0], divisions+1)
y_d = np.linspace(p1[1], p2[1], divisions+1)
z_d = np.linspace(p1[2], p2[2], divisions+1)
for n in xrange(divisions-2):
    p = (x_d[n], y_d[n], z_d[n])
    q = (x_d[n+1], y_d[n+1], z_d[n+1])
    tax.line(p, q, marker=(3, 2, 180+63), **marker_style)

# Text
ax.text(0.36, 0.53, r'$\mathcal{O}_{e\tau}$', fontsize=FONTSIZE,
       rotation=0, verticalalignment='center')
ax.text(0.445, 0.54, r'$\Lambda_d\Longrightarrow$', fontsize=FONTSIZE+5,
       rotation=80, verticalalignment='center')
ax.text(0.68, 0.09, r'$\mathcal{O}_{\mu\tau}$', fontsize=FONTSIZE,
       rotation=0, verticalalignment='center')
ax.text(0.7, 0.14, r'$\Lambda_d\Longrightarrow$', fontsize=FONTSIZE+5,
       rotation=-23, verticalalignment='center')
ax.text(0.34, 0.2, r'$68\%$', fontsize=FONTSIZE, rotation=5)
ax.text(0.34, 0.1, r'$90\%$', fontsize=FONTSIZE, rotation=5)

# Save
outfile = 'fr'
for of in OUTFORMAT:
    o = OUTFOLDER + outfile + '.' + of
    print 'Saving as {0}'.format(o)
    fig.savefig(o, bbox_inches='tight', dpi=150)

# Scenarios
scenarios = [MixingScenario.T12, MixingScenario.T13, MixingScenario.T23]
title = [r'e\mu', r'e\tau', r'\mu\tau']
srcs = [(0, 0, 1), (1, 0, 0), (0, 1, 0), (1, 2, 0)]
colours = {0:PURPLE+[1.], 1:'blue', 2:'green', 3:'red'}
ev_sc = np.linspace(-60, -20, 4E3)

# Figure
fig = plt.figure(figsize=(36, 12))
gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1])
# gs.update(hspace=0.4, wspace=0.15)
for i, sc in enumerate(scenarios):
    ax = fig.add_subplot(gs[i])
    ax.set_aspect('equal')

    # Boundary and Gridlines
    scale = 1
    legend = False
    if i == 1: legend = True
    fig, tax = ternary.figure(ax=ax, scale=scale)
    setupfigure(tax, legend=legend, show_001=True)
    tax.set_title(r'$\mathcal{O}_{'+title[i]+r'}\:operator$',
                  fontsize=FONTSIZE+10, pad=35)

    # Calculate and plot
    for c, s in enumerate(srcs):
        frs = []
        for e in ev_sc:
            u = fr_utils.params_to_BSMu(e, fix_mixing=sc, dim=6, energy=1e6)
            frs.append(fr_utils.u_to_fr(s, np.array(u, dtype=np.complex128)))
        frs = np.vstack(frs)
        tax.scatter(frs, marker='o', s=2, alpha=1, color=colours[c])

    # Contour
    tax.plot(contour_68_upper, linewidth=2.3, color='grey', zorder=0, alpha=0.6)
    tax.plot(contour_68_lower, linewidth=2.3, color='grey', zorder=0, alpha=0.6)
    tax.plot(contour_90_upper, linewidth=2.3, color='darkgrey', zorder=0, alpha=0.6)
    tax.plot(contour_90_lower, linewidth=2.3, color='darkgrey', zorder=0, alpha=0.6)
    ax.text(0.34, 0.2, r'$68\%$', fontsize=FONTSIZE, rotation=5)
    ax.text(0.34, 0.1, r'$90\%$', fontsize=FONTSIZE, rotation=5)

# Save
outfile = 'fr_all'
for of in OUTFORMAT:
    o = OUTFOLDER + outfile + '.' + of
    print 'Saving as {0}'.format(o)
    fig.savefig(o, bbox_inches='tight', dpi=150)
