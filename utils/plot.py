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
from utils.enums import EnergyDependance, Likelihood, ParamTag
from utils.fr import angles_to_u, angles_to_fr

rc('text', usetex=False)
rc('font', **{'family':'serif', 'serif':['Computer Modern'], 'size':18})


def calc_nbins(x):
    n =  (np.max(x) - np.min(x)) / (2 * len(x)**(-1./3) * (np.percentile(x, 75) - np.percentile(x, 25)))
    return np.floor(n)


def most_likely(arr):
    """Return the densest region given a 1D array of data."""
    binning = calc_bins(arr)
    harr = np.histogram(arr, binning)[0]
    return centers(binning)[np.argmax(harr)]


def interval(arr, percentile=68.):
    """Returns the *percentile* shortest interval around the mode."""
    center = most_likely(arr)
    sarr = sorted(arr)
    delta = np.abs(sarr - center)
    curr_low = np.argmin(delta)
    curr_up = curr_low
    npoints = len(sarr)
    while curr_up - curr_low < percentile/100.*npoints:
        if curr_low == 0:
            curr_up += 1
        elif curr_up == npoints-1:
            curr_low -= 1
        elif sarr[curr_up]-sarr[curr_low-1] < sarr[curr_up+1]-sarr[curr_low]:
            curr_low -= 1
        elif sarr[curr_up]-sarr[curr_low-1] > sarr[curr_up+1]-sarr[curr_low]:
            curr_up += 1
        elif (curr_up - curr_low) % 2:
            # they are equal so step half of the time up and down
            curr_low -= 1
        else:
            curr_up += 1
    return sarr[curr_low], center, sarr[curr_up]


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
    t = ''
    mr1, mr2, mr3 = args.measured_ratio
    if args.fix_source_ratio:
        sr1, sr2, sr3 = args.source_ratio
        if args.fix_scale:
            t += 'Source flavour ratio = [{0:.2f}, {1:.2f}, {2:.2f}]\nIC ' \
                    'observed flavour ratio = [{3:.2f}, {4:.2f}, ' \
                    '{5:.2f}]\nDimension = {6}\nScale = {7}'.format(
                        sr1, sr2, sr3, mr1, mr2, mr3, args.dimension,
                        int(args.energy), args.scale
                    )
        else:
            t += 'Source flavour ratio = [{0:.2f}, {1:.2f}, {2:.2f}]\nIC ' \
                    'observed flavour ratio = [{3:.2f}, {4:.2f}, ' \
                    '{5:.2f}]\nDimension = {6}'.format(
                        sr1, sr2, sr3, mr1, mr2, mr3, args.dimension,
                        int(args.energy)
                    )
    else:
        if args.fix_scale:
            t += 'IC observed flavour ratio = [{0:.2f}, {1:.2f}, ' \
                    '{2:.2f}]\nDimension = {3}\nScale = {4}'.format(
                        mr1, mr2, mr3, args.dimension, int(args.energy),
                        args.scale
                    )
	else:
            t += 'IC observed flavour ratio = [{0:.2f}, {1:.2f}, ' \
                    '{2:.2f}]\nDimension = {3}'.format(
                        mr1, mr2, mr3, args.dimension, int(args.energy)
                    )
    if args.likelihood is Likelihood.GAUSSIAN:
        t += '\nSigma = {0:.3f}'.format(args.sigma_ratio)
    if args.energy_dependance is EnergyDependance.SPECTRAL:
        t += '\nSpectral Index = {0}\nBinning = [{1}, {2}] TeV - {3} bins'.format(
            int(args.spectral_index), int(args.binning[0]/1e3),
            int(args.binning[-1]/1e3), len(args.binning)-1
        )
    elif args.energy_dependance is EnergyDependance.MONO:
        t += '\nEnergy = {0} TeV'.format(int(args.energy/1e3))
    return t


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
        nu_index = mcmc_paramset.from_tag(ParamTag.NUISANCE, index=True)
        fr_index = mcmc_paramset.from_tag(ParamTag.MMANGLES, index=True)
        sc_index = mcmc_paramset.from_tag(ParamTag.SCALE, index=True)
        if not args.fix_source_ratio:
            sr_index = mcmc_paramset.from_tag(ParamTag.SRCANGLES, index=True)

        nu_elements = raw[:,nu_index]
        fr_elements = np.array(map(flat_angles_to_u, raw[:,fr_index]))
        sc_elements = raw[:,sc_index]
        if not args.fix_source_ratio:
            sr_elements = np.array(map(angles_to_fr, raw[:,sr_index]))
        if args.fix_source_ratio:
            Tchain = np.column_stack(
                [nu_elements, fr_elements, sc_elements]
            )
        else:
            Tchain = np.column_stack(
                [nu_elements, fr_elements, sc_elements, sr_elements]
            )

        trns_ranges = np.array(ranges)[nu_index,].tolist()
        trns_axes_labels = np.array(axes_labels)[nu_index,].tolist()
        if not args.fix_mixing:
            trns_axes_labels += \
                [r'\mid \tilde{U}_{e1} \mid'    , r'\mid \tilde{U}_{e2} \mid'    , r'\mid \tilde{U}_{e3} \mid'     , \
                 r'\mid \tilde{U}_{\mu1} \mid'  , r'\mid \tilde{U}_{\mu2} \mid'  , r'\mid \tilde{U}_{\mu3} \mid'   , \
                 r'\mid \tilde{U}_{\tau1} \mid' , r'\mid \tilde{U}_{\tau2} \mid' , r'\mid \tilde{U}_{\tau3} \mid']
            trns_ranges += [(0, 1)] * 9
        if not args.fix_scale:
            trns_axes_labels += [np.array(axes_labels)[sc_index].tolist()]
            trns_ranges += [np.array(ranges)[sc_index].tolist()]
        if not args.fix_source_ratio:
            trns_axes_labels += [r'\phi_e', r'\phi_\mu', r'\phi_\tau']
            trns_ranges += [(0, 1)] * 3

        g = plot_Tchain(Tchain, trns_axes_labels, trns_ranges)

        if args.fix_mixing and args.fix_source_ratio:
            mpl.pyplot.figtext(0.4, 0.7, fig_text, fontsize=4)
        else:
            mpl.pyplot.figtext(0.5, 0.7, fig_text, fontsize=15)
        for of in outformat:
            g.export(outfile+'_elements.'+of)
