#! /usr/bin/env python
# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : February 24, 2019

"""
HESE BSM Flavour Figure 2
"""

from __future__ import absolute_import, division

import argparse
from functools import partial

import numpy as np

from utils import fr as fr_utils
from utils import misc as misc_utils
from utils import plot as plot_utils
from utils.plot import PLANCK_SCALE

from matplotlib import pyplot as plt
from matplotlib.patches import Circle
from matplotlib.legend_handler import HandlerPatch
import matplotlib.gridspec as gridspec


DIM = 6
NBINS = 25


class HandlerCircle(HandlerPatch):
    def create_artists(self, legend, orig_handle, xdescent, ydescent, width,
                       height, fontsize, trans):
        r = 10
        x = r + width//2 + 10
        y = height//2 - 3

        # create
        p = Circle(xy=(x, y), radius=r)

        # update with data from original object
        self.update_prop(p, orig_handle, legend)

        # move xy to legend
        p.set_transform(trans)
        return [p]


def cmap_discretize(cmap, N):
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki, key in enumerate(('red','green','blue')):
        cdict[key] = [
            (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in xrange(N+1)
        ]
    # Return colormap object.
    return mpl.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)


def alp(x):
    y = list(x)
    y[-1] = 0.7
    return y


def process_args(args):
    """Process the input args."""
    pass


def parse_args(args=None):
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="HESE BSM Flavour Figure 2",
        formatter_class=misc_utils.SortingHelpFormatter,
    )
    parser.add_argument(
        '--datadir', type=str,
        help='Path to directory containing contour chains'
    )
    if args is None: return parser.parse_args()
    else: return parser.parse_args(args.split())


def main():
    args = parse_args()
    process_args(args)
    misc_utils.print_args(args)

    prefix = ''

    # Load HESE contour.
    contour_infile = args.datadir + '/contour' + prefix + '/contour_REAL.npy'
    contour_angles = np.load(contour_infile)[:,-2:]
    contour_frs = np.array(map(fr_utils.angles_to_fr, contour_angles))

    # Load mc_unitary.
    mcu_basefile = args.datadir + '/mc_unitary' + prefix + '/mc_unitary_SRC_'
    mcu_frs_120 = np.load(mcu_basefile + '1_2_0.npy')
    mcu_frs_100 = np.load(mcu_basefile + '1_0_0.npy')
    mcu_frs_010 = np.load(mcu_basefile + '0_1_0.npy')

    # Load mc_texture.
    mct_basefile = args.datadir + '/mc_texture' + prefix + '/mc_texture_SRC_'
    mct_chains_010_OET = np.load(mct_basefile + '0_1_0_OET.npy')
    mct_chains_100_OUT = np.load(mct_basefile + '1_0_0_OUT.npy')

    # Caculate min/max scales.
    min_scale = 1E100
    max_scale = 1E-100

    for i in (mct_chains_010_OET, mct_chains_100_OUT):
        min_scale = min_scale if min_scale < np.min(i[:,-1]) else np.min(i[:,-1])
        max_scale = max_scale if max_scale > np.max(i[:,-1]) else np.max(i[:,-1])
    min_scale -= np.log10(PLANCK_SCALE[DIM])
    max_scale -= np.log10(PLANCK_SCALE[DIM])

    cmap_green = cmap_discretize(mpl.colors.LinearSegmentedColormap.from_list(
        "", ["lime", "gold", "darkorange"]
    ), 10)
    cmap_blue = cmap_discretize(mpl.colors.LinearSegmentedColormap.from_list(
        "", ["blue", "fuchsia", "darkmagenta"]
    ), 10)

    norm = mpl.colors.Normalize(vmin=min_scale, vmax=max_scale)
    color_010 = map(alp, map(cmap_green, map(
        norm, mct_chains_010_OET[:,-1]-np.log10(PLANCK_SCALE[DIM])
    )))
    color_100 = map(alp, map(cmap_blue, map(
        norm, mct_chains_100_OUT[:,-1]-np.log10(PLANCK_SCALE[DIM])
    )))

    fontsize = 23
    ax_labels = [r'$f_{e}^{\oplus}$', r'$f_{\mu}^{\oplus}$', r'$f_{\tau}^{\oplus}$']

    # Figure
    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[40, 1])
    gs.update(hspace=0.3, wspace=0.15)

    # Axis
    ax = fig.add_subplot(gs[0])
    tax = plot_utils.get_tax(ax, scale=NBINS, ax_labels=ax_labels)

    # Plot HESE contour
    coverages = {68: 'grey', 90: 'darkgrey'}
    for cov in coverages.iterkeys():
        plot_utils.flavour_contour(
            frs = contour_frs,
            ax = ax,
            nbins = NBINS,
            coverage = cov,
            linewidth = 2.3,
            color = coverages[cov],
            alpha=0.6,
            zorder=0
        )
    ax.text(0.34*NBINS, 0.143*NBINS, r'$68\%$', fontsize=fontsize, rotation=3)
    ax.text(0.34*NBINS, 0.038*NBINS, r'$90\%$', fontsize=fontsize, rotation=0)

    # Plot unitary contours
    mcu_kwargs = {mcu_frs_120:{'color':'red',   'alpha':0.2, 'zorder':3},
                  mcu_frs_100:{'color':'black', 'alpha':0.1, 'zorder':2},
                  mcu_frs_010:{'color':'black', 'alpha':0.1, 'zorder':2}}
    for frs, kwargs in mcu_kwargs.itervals():
        plot_utils.flavour_contour(
            frs = frs,
            ax = ax,
            fill = True,
            nbins = NBINS,
            coverage = 100,
            linewidth = 1,
            **kwargs
        )

    # Plot BSM points
    tax.scatter(
        mct_chains_010_OET[:,:-1]*NBINS, marker='o', s=0.03, color=color_010,
        zorder=5
    )
    tax.scatter(
        mct_chains_100_OUT[:,:-1]*NBINS, marker='o', s=0.03, color=color_100,
        zorder=5
    )

    # Legend
    legend_elements = []
    legend_elements.append(
        Circle((0., 0.), 0.1, facecolor='lime', alpha=0.7, edgecolor='k',
               linewidth=2., label=r'$\left (0:1:0\right )\:w/\:{\rm New\:Physics}$')
    )
    legend_elements.append(
        Circle((0., 0.), 0.1, facecolor='blue', alpha=0.7, edgecolor='k',
               linewidth=2., label=r'$\left (1:0:0\right )\:w/\:{\rm New\:Physics}$')
    )
    legend_elements.append(
        Circle((0., 0.), 0.1, facecolor='red', alpha=0.7, edgecolor='k',
               linewidth=2., label=r'$\left (1:2:0\right )$')
    )
    legend_elements.append(
        Circle((0., 0.), 0.1, facecolor='grey', alpha=0.7, edgecolor='k',
               linewidth=2., label=r'$\left (0:1:0\right ) + \left (1:0:0\right )$')
    )
    legend = plt.legend(handles=legend_elements, loc=(0.65, 0.8),
                        title='Source composition',
                        fontsize=fontsize,
                        handler_map={Circle: HandlerCircle()})
    plt.setp(legend.get_title(), fontsize=fontsize)
    legend.get_frame().set_linestyle('-')

    # Colorbar
    gs00 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1])
    ax0 = fig.add_subplot(gs00[0])
    cb = mpl.colorbar.ColorbarBase(
        ax0, cmap=cmap_010, norm=norm, orientation='horizontal'
    )
    cb.ax.tick_params(labelsize=fontsize-5)
    ax0.text(
        0.5, 2, r'$\mathcal{O}_{e\tau}\:texture$', fontsize=fontsize,
        rotation=0, verticalalignment='center', horizontalalignment='center'
    )


    ax1 = fig.add_subplot(gs00[1])
    cb = mpl.colorbar.ColorbarBase(ax1, cmap=cmap_100, norm=norm, orientation='horizontal')
    cb.ax.tick_params(labelsize=fontsize-5)
    ax1.text(0.5, 2, r'$\mathcal{O}_{\mu\tau}\:texture$', fontsize=fontsize,
             rotation=0, verticalalignment='center', horizontalalignment='center')

    fig.text(0.5, 0.038, r'${\rm New\:Physics\:Scale}\:[\:{\rm log}_{10} (\Lambda_6\:/\:{\rm M}^{2}_{\rm Planck})\: ]$',
             fontsize=fontsize+5, horizontalalignment='center')

    outformat = ['png']
    outfile = args.datadir[:5]+args.datadir[5:].replace('data', 'plots')
    outfile += '/fig2' + prefix
    make_dir(outfile)
    for of in outformat:
        print 'Saving plot as {0}'.format(outfile + '.' + of)
        fig.savefig(outfile + '.' + of, bbox_inches='tight', dpi=150)

    print "DONE!"


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
