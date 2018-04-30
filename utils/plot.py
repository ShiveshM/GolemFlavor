# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : March 19, 2018

"""
Plotting functions for the BSM flavour ratio analysis
"""

from __future__ import absolute_import, division

import os

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import rc
from matplotlib import pyplot as plt
from matplotlib.offsetbox import AnchoredText

from getdist import plots, mcsamples

from utils import misc as misc_utils
from utils.enums import EnergyDependance, Likelihood, ParamTag, StatCateg
from utils.fr import angles_to_u, angles_to_fr

rc('text', usetex=False)
rc('font', **{'family':'serif', 'serif':['Computer Modern'], 'size':18})


def centers(x):
    return (x[:-1]+x[1:])*0.5


def get_units(dimension):
    if dimension == 3: return r' / GeV'
    if dimension == 4: return r''
    if dimension == 5: return r' / GeV^{-1}'
    if dimension == 6: return r' / GeV^{-2}'
    if dimension == 7: return r' / GeV^{-3}'
    if dimension == 8: return r' / GeV^{-4}'


def calc_nbins(x):
    n =  (np.max(x) - np.min(x)) / (2 * len(x)**(-1./3) * (np.percentile(x, 75) - np.percentile(x, 25)))
    return np.floor(n)


def calc_bins(x):
    nbins = calc_nbins(x)
    return np.linspace(np.min(x), np.max(x)+2, num=nbins+1)


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
        t += 'Source flavour ratio = [{0:.2f}, {1:.2f}, {2:.2f}]\nIC ' \
                'observed flavour ratio = [{3:.2f}, {4:.2f}, ' \
                '{5:.2f}]\nDimension = {6}'.format(
                    sr1, sr2, sr3, mr1, mr2, mr3, args.dimension,
                    int(args.energy)
                )
    else:
        t += 'IC observed flavour ratio = [{0:.2f}, {1:.2f}, ' \
                '{2:.2f}]\nDimension = {3}'.format(
                    mr1, mr2, mr3, args.dimension, int(args.energy)
                )
    if args.fix_scale:
        t += 'Scale = {0}'.format(args.scale)
    if args.likelihood is Likelihood.GAUSSIAN:
        t += '\nSigma = {0:.3f}'.format(args.sigma_ratio)
    if args.energy_dependance is EnergyDependance.SPECTRAL:
        if not args.fold_index:
            t += '\nSpectral Index = {0}'.format(int(args.spectral_index))
        t += '\nBinning = [{0}, {1}] TeV - {2} bins'.format(
            int(args.binning[0]/1e3), int(args.binning[-1]/1e3), len(args.binning)-1
        )
    elif args.energy_dependance is EnergyDependance.MONO:
        t += '\nEnergy = {0} TeV'.format(int(args.energy/1e3))
    return t


def chainer_plot(infile, outfile, outformat, args, llh_paramset):
    """Make the triangle plot."""
    if not args.plot_angles and not args.plot_elements:
        return

    raw = np.load(infile)
    print 'raw.shape', raw.shape

    misc_utils.make_dir(outfile)
    fig_text = gen_figtext(args)

    axes_labels = llh_paramset.labels
    ranges = llh_paramset.ranges

    if args.plot_angles:
        print "Making triangle plots"
        Tchain = raw
        g = plot_Tchain(Tchain, axes_labels, ranges)

        mpl.pyplot.figtext(0.5, 0.7, fig_text, fontsize=15)

        for i_ax_1, ax_1 in enumerate(g.subplots):
            for i_ax_2, ax_2 in enumerate(ax_1):
                if i_ax_1 == i_ax_2:
                    itv = interval(Tchain[:,i_ax_1], percentile=90.)
                    for l in itv:
                        ax_2.axvline(l, color='gray', ls='--')
                        ax_2.set_title(r'${0:.2f}_{{{1:.2f}}}^{{+{2:.2f}}}$'.format(
                            itv[1], itv[0]-itv[1], itv[2]-itv[1]
                        ), fontsize=10)

        # if not args.fix_mixing:
        #     sc_index = llh_paramset.from_tag(ParamTag.SCALE, index=True)
        #     itv = interval(Tchain[:,sc_index], percentile=90.)
        #     mpl.pyplot.figtext(
        #         0.5, 0.3, 'Scale 90% Interval = [1E{0}, 1E{1}], Center = '
        #         '1E{2}'.format(itv[0], itv[2], itv[1])
        #     )

        for of in outformat:
            g.export(outfile+'_angles.'+of)

    if args.plot_elements:
        print "Making triangle plots"
        if args.fix_mixing_almost:
            raise NotImplementedError
        nu_index = llh_paramset.from_tag(ParamTag.NUISANCE, index=True)
        fr_index = llh_paramset.from_tag(ParamTag.MMANGLES, index=True)
        sc_index = llh_paramset.from_tag(ParamTag.SCALE, index=True)
        if not args.fix_source_ratio:
            sr_index = llh_paramset.from_tag(ParamTag.SRCANGLES, index=True)

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

        mpl.pyplot.figtext(0.5, 0.7, fig_text, fontsize=15)
        for of in outformat:
            g.export(outfile+'_elements.'+of)


def plot_statistic(data, outfile, outformat, args, scale_param, label=None):
    """Make MultiNest factor or LLH value plot."""
    print "Making Statistic plot"
    fig_text = gen_figtext(args)
    if label is not None: fig_text += '\n' + label

    print 'data', data
    print 'data.shape', data.shape
    scales, statistic = data.T
    if args.stat_method is StatCateg.BAYESIAN:
        min_idx = np.argmin(scales)
        null = statistic[min_idx]
        reduced_ev = -(statistic - null)
    elif args.stat_method is StatCateg.FREQUENTIST:
        min_idx = np.argmin(scales)
        null = statistic[min_idx]
        reduced_ev = -2*(statistic - null)

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)

    ax.set_xlim(scale_param.ranges)
    ax.set_xlabel('$'+scale_param.tex+'$')
    if args.stat_method is StatCateg.BAYESIAN:
        ax.set_ylabel(r'Bayes Factor')
    elif args.stat_method is StatCateg.FREQUENTIST:
        ax.set_ylabel(r'$-2\Delta {\rm LLH}$')

    ax.plot(scales, reduced_ev)

    for ymaj in ax.yaxis.get_majorticklocs():
        ax.axhline(y=ymaj, ls=':', color='gray', alpha=0.3, linewidth=1)
    for xmaj in ax.xaxis.get_majorticklocs():
        ax.axvline(x=xmaj, ls=':', color='gray', alpha=0.3, linewidth=1)

    at = AnchoredText(
        '\n'+fig_text, prop=dict(size=7), frameon=True, loc=2
    )
    at.patch.set_boxstyle("round,pad=0.1,rounding_size=0.5")
    ax.add_artist(at)
    misc_utils.make_dir(outfile)
    for of in outformat:
        fig.savefig(outfile+'.'+of, bbox_inches='tight', dpi=150)


def myround(x, base=5, up=False, down=False):
    if up == down and up is True: assert 0
    if up: return int(base * np.round(float(x)/base-0.5))
    elif down: return int(base * np.round(float(x)/base+0.5))
    else: int(base * np.round(float(x)/base))


def plot_BSM_angles_limit(dirname, outfile, outformat, args, bayesian):
    """Make BSM angles vs scale limit plot."""
    if not args.plot_angles_limit: return
    print "Making BSM angles limit plot."""
    fig_text = gen_figtext(args)
    xticks = [r'$\mathcal{O}_{12}$', r'$\mathcal{O}_{13}$', r'$\mathcal{O}_{23}$']

    raw = []
    for root, dirs, filenames in os.walk(dirname):
        for fn in filenames:
            if fn[-4:] == '.npy':
                raw.append(np.load(os.path.join(root, fn)))
    raw = np.vstack(raw)
    raw = raw[np.argsort(raw[:,0])]
    print 'raw', raw
    print 'raw.shape', raw.shape
    sc_ranges = (
        myround(np.min(raw[0][:,0]), up=True),
        myround(np.max(raw[0][:,0]), down=True)
    )

    proc = []
    if bayesian:
        for idx, theta in enumerate(raw):
            scale, evidences = theta.T
            null = evidences[0]
            reduced_ev = -(evidences - null)
            al = scale[reduced_ev > np.log(10**(1/2.))]
            proc.append((idx+1, al[0]))
    else:
        for idx, theta in enumerate(raw):
            scale, llh = theta.T
            delta_llh = -2 * (llh - np.max(llh))
            # 90% CL for 1 dof
            al = scale[delta_llh > 2.71]
            proc.append((idx+1, al[0]))

    limits = np.array(proc)
    print 'limits',  limits

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)

    ax.set_xticklabels(['']+xticks+[''])
    ax.set_xlim(0, len(xticks)+1)
    ax.set_ylim(sc_ranges[0], sc_ranges[-1])
    ax.set_xlabel(r'BSM angle')
    ylabel = r'${\rm log}_{10} \Lambda' + get_units(args.dimension) + r'$'
    ax.set_ylabel(ylabel)

    for l in limits:
        line = plt.Line2D(
            (l[0]-0.1, l[0]+0.1), (l[1], l[1]), lw=3, color='r'
        )
        ax.add_line(line)
        # ax.arrow(
        #     l[0], l[1], 0, -1.5, head_width=0.05, head_length=0.2, fc='r',
        #     ec='r', lw=2
        # )
        ax.annotate(
            s='', xy=l, xytext=(l[0], l[1]+1.5),
            arrowprops={'arrowstyle': '<-', 'lw': 1.5, 'color':'r'}
        )

    for ymaj in ax.yaxis.get_majorticklocs():
        ax.axhline(y=ymaj, ls='-', color='gray', alpha=0.4, linewidth=1)
    for xmaj in ax.xaxis.get_majorticklocs():
        ax.axvline(x=xmaj, ls='-', color='gray', alpha=0.4, linewidth=1)

    for of in outformat:
        fig.savefig(outfile+'.'+of, bbox_inches='tight', dpi=150)

