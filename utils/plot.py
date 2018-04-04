# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : March 19, 2018

"""
Plotting functions for the BSM flavour ratio analysis
"""

from __future__ import absolute_import, division

import argparse

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import rc

import getdist
from getdist import plots
from getdist import mcsamples

from utils import misc as misc_utils
from utils.enums import ParamTag
from utils.fr import angles_to_u, angles_to_fr

rc('text', usetex=False)
rc('font', **{'family':'serif', 'serif':['Computer Modern'], 'size':18})


def plot_argparse(parser):
    """Arguments for plotting."""
    parser.add_argument(
        '--plot-angles', type=misc_utils.parse_bool, default='True',
        help='Plot MCMC triangle in the angles space'
    )
    parser.add_argument(
        '--plot-elements', type=misc_utils.parse_bool, default='False',
        help='Plot MCMC triangle in the mixing elements space'
    )


def flat_angles_to_u(x):
    """Convert from angles to mixing elements."""
    return abs(angles_to_u(x)).astype(np.float32).flatten().tolist()


def plot_Tchain(Tchain, axes_labels, ranges):
    """Plot the Tchain using getdist."""
    Tsample = mcsamples.MCSamples(
        samples=Tchain, labels=axes_labels, ranges=ranges
    )

    Tsample.updateSettings({'contours': [0.90, 0.99]})
    Tsample.num_bins_2D=500
    Tsample.fine_bins_2D=500
    Tsample.smooth_scale_2D=0.03

    g = plots.getSubplotPlotter()
    g.settings.num_plot_contours = 2
    g.settings.axes_fontsize = 10
    g.settings.figure_legend_frame = False
    g.triangle_plot(
        [Tsample], filled=True,
    )
    return g


def gen_figtext(args):
    """Generate the figure text."""
    mr1, mr2, mr3 = args.measured_ratio
    if args.fix_source_ratio:
        sr1, sr2, sr3 = args.source_ratio
        if args.fix_scale:
            return 'Source flavour ratio = [{0:.2f}, {1:.2f}, {2:.2f}]\nIC ' \
                    'observed flavour ratio = [{3:.2f}, {4:.2f}, ' \
                    '{5:.2f}]\nSigma = {6:.3f}\nDimension = {7}\nEnergy = ' \
                    '{8} GeV\nScale = {9}'.format(
                        sr1, sr2, sr3, mr1, mr2, mr3, args.sigma_ratio,
                        args.dimension, int(args.energy), args.scale
                    )
        else:
            return 'Source flavour ratio = [{0:.2f}, {1:.2f}, {2:.2f}]\nIC ' \
                    'observed flavour ratio = [{3:.2f}, {4:.2f}, ' \
                    '{5:.2f}]\nSigma = {6:.3f}\nDimension = {7}\nEnergy = {8} ' \
                    'GeV'.format(
                        sr1, sr2, sr3, mr1, mr2, mr3, args.sigma_ratio,
                        args.dimension, int(args.energy)
                    )
    else:
        if args.fix_scale:
            return 'IC observed flavour ratio = [{0:.2f}, {1:.2f}, ' \
                    '{2:.2f}]\nSigma = {3:.3f}\nDimension = {4}\nEnergy = {5} ' \
                    'GeV\nScale = {6}'.format(
                        mr1, mr2, mr3, args.sigma_ratio, args.dimension,
                        int(args.energy), args.scale
                    )
	else:
            return 'IC observed flavour ratio = [{0:.2f}, {1:.2f}, ' \
                    '{2:.2f}]\nSigma = {3:.3f}\nDimension = {4}\nEnergy = {5} ' \
                    'GeV'.format(
                        mr1, mr2, mr3, args.sigma_ratio, args.dimension,
                        int(args.energy)
                    )

def chainer_plot(infile, outfile, outformat, args, mcmc_paramset):
    """Make the triangle plot."""
    if not args.plot_angles and not args.plot_elements:
        return

    raw = np.load(infile)
    print 'raw.shape', raw.shape

    misc_utils.make_dir(outfile)
    fig_text = gen_figtext(args)

    axes_labels = mcmc_paramset.labels
    ranges = mcmc_paramset.ranges

    if args.plot_angles:
        Tchain = raw
        g = plot_Tchain(Tchain, axes_labels, ranges)

        if args.fix_mixing and args.fix_source_ratio:
            mpl.pyplot.figtext(0.4, 0.7, fig_text, fontsize=4)
        else:
            mpl.pyplot.figtext(0.5, 0.7, fig_text, fontsize=15)
        for of in outformat:
            g.export(outfile+'_angles.'+of)

    if args.plot_elements:
        nu_index = mcmc_paramset.idx_from_tag(ParamTag.NUISANCE)
        fr_index = mcmc_paramset.idx_from_tag(ParamTag.MMANGLES)
        sc_index = mcmc_paramset.idx_from_tag(ParamTag.SCALE)
        sr_index = mcmc_paramset.idx_from_tag(ParamTag.SRCANGLES)

        nu_elements = raw[:,nu_index]
        fr_elements = np.array(map(angles_to_fr, raw[:,fr_index]))
        sc_elements = raw[:,sc_index]
        sr_elements = np.array(map(flat_angles_to_u, raw[:,sr_index]))
        Tchain = np.column_stack(
            [nu_elements, fr_elements, sc_elements, sr_elements]
        )

        trns_ranges = ranges[nu_index,]
        trns_axes_labels = axes_labels[nu_index,]
        if not args.fix_mixing:
            trns_axes_labels += \
                [r'\mid \tilde{U}_{e1} \mid'    , r'\mid \tilde{U}_{e2} \mid'    , r'\mid \tilde{U}_{e3} \mid'     , \
                 r'\mid \tilde{U}_{\mu1} \mid'  , r'\mid \tilde{U}_{\mu2} \mid'  , r'\mid \tilde{U}_{\mu3} \mid'   , \
                 r'\mid \tilde{U}_{\tau1} \mid' , r'\mid \tilde{U}_{\tau2} \mid' , r'\mid \tilde{U}_{\tau3} \mid']
            trns_ranges += [(0, 1)] * 9
        if not args.fix_scale:
            trns_axes_labels += [axes_labels[sc_index]]
            trns_ranges += [ranges[sc_index]]
        if not args.fix_source_ratio:
            trns_axes_labels += [r'\phi_e', r'\phi_\mu', r'\phi_\tau']
            trns_ranges += [(0, 1)] * 3

        g = plot_Tchain(Tchain, trns_axes_labels, trns_ranges)

        if args.fix_mixing and args.fix_source_ratio:
            mpl.pyplot.figtext(0.4, 0.7, fig_text, fontsize=4)
        else:
            mpl.pyplot.figtext(0.5, 0.7, fig_text, fontsize=15)
        for of in outformat:
            g.export(outfile+'_angles.'+of)
