# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : March 19, 2018

"""
Plotting functions for the BSM flavour ratio analysis
"""

from __future__ import absolute_import, division

import os
import socket
from copy import deepcopy

import numpy as np
from scipy import interpolate

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import rc
from matplotlib import pyplot as plt
from matplotlib.offsetbox import AnchoredText
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from getdist import plots, mcsamples

from utils import misc as misc_utils
from utils.enums import DataType, EnergyDependance
from utils.enums import Likelihood, ParamTag, StatCateg
from utils.fr import angles_to_u, angles_to_fr

plt.style.use(os.environ['GOLEMSOURCEPATH']+'/GolemFit/scripts/paper/paper.mplstyle')
if 'submitter' in socket.gethostname():
    rc('text', usetex=False)


def centers(x):
    return (x[:-1]+x[1:])*0.5


def get_units(dimension):
    if dimension == 3: return r' / \:{\rm GeV}'
    if dimension == 4: return r''
    if dimension == 5: return r' / \:{rm GeV}^{-1}'
    if dimension == 6: return r' / \:{rm GeV}^{-2}'
    if dimension == 7: return r' / \:{rm GeV}^{-3}'
    if dimension == 8: return r' / \:{rm GeV}^{-4}'


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
    mr1, mr2, mr3 = misc_utils.solve_ratio(args.measured_ratio)
    if args.fix_source_ratio:
        sr1, sr2, sr3 = misc_utils.solve_ratio(args.source_ratio)
        t += 'Source flavour ratio = [{0}, {1}, {2}]'.format(sr1, sr2, sr3)
        if args.data in [DataType.ASIMOV, DataType.REALISATION]:
            t += '\nIC observed flavour ratio = [{0}, {1}, {2}]'.format(
                mr1, mr2, mr3
            )
        t += '\nDimension = {0}'.format(args.dimension)
    else:
        if args.data in [DataType.ASIMOV, DataType.REALISATION]:
            t += 'IC observed flavour ratio = [{0}, {1}, ' \
                    '{2}]\nDimension = {3}'.format(
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

        if args.data is DataType.REAL:
            fig.text(0.8, 0.7, 'IceCube Preliminary', color='red', fontsize=15,
                     ha='center', va='center')
        elif args.data in [DataType.ASIMOV, DataType.REALISATION]:
            fig.text(0.8, 0.7, 'IceCube Simulation', color='red', fontsize=15,
                     ha='center', va='center')

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

        if args.data is DataType.REAL:
            fig.text(0.8, 0.7, 'IceCube Preliminary', color='red', fontsize=15,
                     ha='center', va='center')
        elif args.data in [DataType.ASIMOV, DataType.REALISATION]:
            fig.text(0.8, 0.7, 'IceCube Simulation', color='red', fontsize=15,
                     ha='center', va='center')

        mpl.pyplot.figtext(0.5, 0.7, fig_text, fontsize=15)
        for of in outformat:
            g.export(outfile+'_elements.'+of)


def myround(x, base=2, up=False, down=False):
    if up == down and up is True: assert 0
    if up: return int(base * np.round(float(x)/base-0.5))
    elif down: return int(base * np.round(float(x)/base+0.5))
    else: int(base * np.round(float(x)/base))


def plot_statistic(data, outfile, outformat, args, scale_param, label=None):
    """Make MultiNest factor or LLH value plot."""
    print 'Making Statistic plot'
    fig_text = gen_figtext(args)
    if label is not None: fig_text += '\n' + label

    print 'data', data
    print 'data.shape', data.shape
    scales, statistic = data.T
    tck, u = interpolate.splprep([scales, statistic], s=0)
    scales, statistic = interpolate.splev(np.linspace(0, 1, 1000), tck)
    print 'scales', scales
    print 'statistic', statistic

    min_idx = np.argmin(scales)
    null = statistic[min_idx]
    if args.stat_method is StatCateg.BAYESIAN:
        reduced_ev = -(statistic - null)
    elif args.stat_method is StatCateg.FREQUENTIST:
        reduced_ev = -2*(statistic - null)

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)

    ax.set_xlim(np.log10(args.scale_region))
    ax.set_xlabel(r'${\mathrm {log}}_{10} \left (\Lambda^{-1}' + \
                  get_units(args.dimension) +r'\right )$', fontsize=16)
    if args.stat_method is StatCateg.BAYESIAN:
        ax.set_ylabel(r'log(Bayes Factor)')
    elif args.stat_method is StatCateg.FREQUENTIST:
        ax.set_ylabel(r'$-2\Delta {\rm LLH}$')

    ax.plot(scales, reduced_ev)

    ax.axhline(y=np.log(10**(3/2.)), color='red', alpha=1., linewidth=1.3)

    for ymaj in ax.yaxis.get_majorticklocs():
        ax.axhline(y=ymaj, ls=':', color='gray', alpha=0.3, linewidth=1)
    for xmaj in ax.xaxis.get_majorticklocs():
        ax.axvline(x=xmaj, ls=':', color='gray', alpha=0.3, linewidth=1)

    if args.data is DataType.REAL:
        fig.text(0.8, 0.14, 'IceCube Preliminary', color='red', fontsize=9,
                 ha='center', va='center')
    elif args.data in [DataType.ASIMOV, DataType.REALISATION]:
        fig.text(0.8, 0.14, 'IceCube Simulation', color='red', fontsize=9,
                 ha='center', va='center')

    at = AnchoredText(
        fig_text, prop=dict(size=10), frameon=True, loc=4
    )
    at.patch.set_boxstyle("round,pad=0.1,rounding_size=0.5")
    ax.add_artist(at)
    misc_utils.make_dir(outfile)
    for of in outformat:
        print 'Saving as {0}'.format(outfile+'.'+of)
        fig.savefig(outfile+'.'+of, bbox_inches='tight', dpi=150)


def plot_sens_full(data, outfile, outformat, args):
    print 'Making FULL sensitivity plot'

    colour = {0:'red', 1:'blue', 2:'green', 3:'purple', 4:'orange', 5:'black'}
    xticks = ['{0}'.format(x) for x in range(np.min(args.dimensions),
                                             np.max(args.dimensions)+1)]
    yranges = [np.inf, -np.inf]
    legend_handles = []

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)

    ax.set_xlim(args.dimensions[0]-1, args.dimensions[-1]+1)
    ax.set_xticklabels([''] + xticks + [''])
    ax.set_xlabel(r'BSM operator dimension ' + r'$d$')
    ax.set_ylabel(r'${\rm log}_{10} \left (\Lambda^{-1} / GeV^{-d+4} \right )$')

    argsc = deepcopy(args)
    for idim in xrange(len(data)):
        dim = args.dimensions[idim]
        print 'dim', dim
        argsc.dimension = dim
        for isrc in xrange(len(data[idim])):
            src = args.source_ratios[isrc]
            argsc.source_ratio = src

            # fig_text = gen_figtext(argsc)
            scales, statistic = data[idim][isrc].T
            min_idx = np.argmin(scales)
            null = statistic[min_idx]
            if args.stat_method is StatCateg.BAYESIAN:
                reduced_ev = -(statistic - null)
                al = scales[reduced_ev > np.log(10**(3/2.))] # Strong degree of belief
                # al = scales[reduced_ev > 0.4] # Testing
            elif args.stat_method is StatCateg.FREQUENTIST:
                reduced_ev = -2*(statistic - null)
                al = scales[reduced_ev > 2.71] # 90% CL for 1 DOF via Wilks
            if len(al) == 0:
                print 'No points for DIM {0} FRS {1} NULL {2}!'.format(
                    dim, misc_utils.solve_ratio(src), null
                )
                print 'Reduced EV {0}'.format(reduced_ev)
                continue
            lim = al[0]
            label = '[{0}, {1}, {2}]'.format(*misc_utils.solve_ratio(src))
            if lim < yranges[0]: yranges[0] = lim
            if lim > yranges[1]: yranges[1] = lim+4
            line = plt.Line2D(
                (dim-0.1, dim+0.1), (lim, lim), lw=3, color=colour[isrc], label=label
            )
            ax.add_line(line)
            if idim == 0: legend_handles.append(line)
            x_offset = isrc*0.05 - 0.05
            ax.annotate(
                s='', xy=(dim+x_offset, lim), xytext=(dim+x_offset, lim+3),
                arrowprops={'arrowstyle': '<-', 'lw': 1.2, 'color':colour[isrc]}
            )

    try:
        yranges = (myround(yranges[0], up=True), myround(yranges[1], down=True))
        ax.set_ylim(yranges)
    except: pass

    ax.legend(handles=legend_handles, prop=dict(size=8), loc='upper right')
    for ymaj in ax.yaxis.get_majorticklocs():
        ax.axhline(y=ymaj, ls=':', color='gray', alpha=0.4, linewidth=1)
    for xmaj in ax.xaxis.get_majorticklocs():
        ax.axvline(x=xmaj, ls=':', color='gray', alpha=0.4, linewidth=1)

    misc_utils.make_dir(outfile)
    for of in outformat:
        print 'Saving plot as {0}'.format(outfile+'.'+of)
        fig.savefig(outfile+'.'+of, bbox_inches='tight', dpi=150)


def plot_sens_fixed_angle_pretty(data, outfile, outformat, args):
    print 'Making FIXED_ANGLE_PRETTY sensitivity plot'
    argsc = deepcopy(args)
    dims = len(data)
    LV_atmo_90pc_limits = {
        3: (2E-24, 1E-1),
        4: (2.7E-28, 3.16E-25),
        5: (1.5E-32, 1.12E-27),
        6: (9.1E-37, 2.82E-30),
        7: (3.6E-41, 1.77E-32),
        8: (1.4E-45, 1.00E-34)
    }

    show_data = True

    en = np.log10([1E4, 1E7])
    bote = {
        3: (-21-(en[0]+en[0]*0), -21-(en[1]+en[1]*0)),
        4: (-21-(en[0]+en[0]*1), -21-(en[1]+en[1]*1)),
        5: (-21-(en[0]+en[0]*2), -21-(en[1]+en[1]*2)),
        6: (-21-(en[0]+en[0]*3), -21-(en[1]+en[1]*3)),
        7: (-21-(en[0]+en[0]*4), -21-(en[1]+en[1]*4)),
        8: (-21-(en[0]+en[0]*5), -21-(en[1]+en[1]*5))
    }

    colour = {0:'red', 1:'blue', 2:'green'}
    rgb_co = {0:[1,0,0], 1:[0,0,1], 2:[0.0, 0.5019607843137255, 0.0]}
    yticks = [r'$\mathcal{O}_{e\mu}$', r'$\mathcal{O}_{e\tau}$', r'$\mathcal{O}_{\mu\tau}$']
    xlims = (-60, -20)

    fig = plt.figure(figsize=(8, 6))
    gs = gridspec.GridSpec(dims, 1)
    gs.update(hspace=0.15)

    first_ax = None
    legend_log = []
    legend_elements = []

    for idim in xrange(len(data)):
        dim = args.dimensions[idim]
        print '== dim', dim
        argsc.dimension = dim
        gs0 = gridspec.GridSpecFromSubplotSpec(
            len(yticks), 1, subplot_spec=gs[idim], hspace=0
        )

        for ian in xrange(len(yticks)):
            print '=== an', ian
            ax = fig.add_subplot(gs0[ian])
            if first_ax is None: first_ax = ax
            ax.set_xlim(xlims)
            ylim = (0.5, 1.5)
            ax.set_ylim(ylim)
            # ax.set_yticks([ylim[0], 1., ylim[1]])
            # ax.set_yticklabels([''] + [yticks[ian]] + [''], fontsize=13)
            # ax.yaxis.get_major_ticks()[0].set_visible(False)
            # ax.yaxis.get_major_ticks()[2].set_visible(False)
            ax.set_yticks([1.])
            ax.set_yticklabels([yticks[ian]], fontsize=13)
            ax.yaxis.tick_right()
            for xmaj in ax.xaxis.get_majorticklocs():
                ax.axvline(x=xmaj, ls=':', color='gray', alpha=0.2, linewidth=1)
            ax.get_xaxis().set_visible(False)
            linestyle = (5, (10, 10))
            if ian == 0:
                ax.spines['bottom'].set_alpha(0.6)
            elif ian == 1:
                ax.text(
                    -0.04, ylim[0], '$d = {0}$'.format(dim), fontsize=16,
                    rotation='90', verticalalignment='center',
                    transform=ax.transAxes
                )
                dim_label_flag = False
                ax.spines['top'].set_alpha(0.6)
                ax.spines['bottom'].set_alpha(0.6)
            elif ian == 2:
                ax.spines['top'].set_alpha(0.6)

            for isrc in xrange(len(data[idim])):
                src = args.source_ratios[isrc]
                print '== src', src
                argsc.source_ratio = src

                if show_data:
                    alpha = 0.03
                    col = 'black'
                else:
                    alpha = 0.07
                    col = 'blue'
                ax.add_patch(patches.Rectangle(
                    (bote[dim][1], ylim[0]), bote[dim][0]-bote[dim][1], np.diff(ylim),
                    fill=True, facecolor=col, alpha=alpha, linewidth=0
                ))

                scales, statistic = data[idim][isrc][ian].T
                tck, u = interpolate.splprep([scales, statistic], s=0)
                scales, statistic = interpolate.splev(np.linspace(0, 1, 1000), tck)
                min_idx = np.argmin(scales)
                null = statistic[min_idx]
                if args.stat_method is StatCateg.BAYESIAN:
                    reduced_ev = -(statistic - null)
                    al = scales[reduced_ev > np.log(10**(3/2.))] # Strong degree of belief
                elif args.stat_method is StatCateg.FREQUENTIST:
                    reduced_ev = -2*(statistic - null)
                    al = scales[reduced_ev > 2.71] # 90% CL for 1 DOF via Wilks
                if len(al) == 0:
                    print 'No points for DIM {0} FRS {1}!'.format(dim, src)
                    continue
                if reduced_ev[-1] < np.log(10**(3/2.)) - 0.1:
                    print 'Peaked contour does not exclude large scales! For ' \
                        'DIM {0} FRS{1}!'.format(dim, src)
                    continue
                lim = al[0]
                print 'limit = {0}'.format(lim)

                if show_data:
                    ax.axvline(x=lim, color=colour[isrc], alpha=1., linewidth=1.5)
                    ax.add_patch(patches.Rectangle(
                        (lim, ylim[0]), 100, np.diff(ylim), fill=True, facecolor=colour[isrc],
                        alpha=0.3, linewidth=0
                    ))

                    if isrc not in legend_log:
                        legend_log.append(isrc)
                        label = '{0} : {1} : {2} at source'.format(*misc_utils.solve_ratio(src))
                        legend_elements.append(
                            Patch(facecolor=rgb_co[isrc]+[0.3],
                                  edgecolor=rgb_co[isrc]+[1], label=label)
                        )

            if ian == 2:
                LV_lim = np.log10(LV_atmo_90pc_limits[dim])
                ax.add_patch(patches.Rectangle(
                    (LV_lim[1], ylim[0]), LV_lim[0]-LV_lim[1], np.diff(ylim),
                    fill=False, hatch='\\\\'
                ))

    ax.get_xaxis().set_visible(True)
    ax.set_xlabel(r'${\rm New\:Physics\:Scale}\:[\:{\rm log}_{10} (\Lambda^{-1}\:/\:{\rm GeV}^{-d+4})\: ]$',
                 fontsize=19)
    ax.tick_params(axis='x', labelsize=16)

    legend_elements.append(
        Patch(facecolor=col, alpha=alpha+0.1, edgecolor='k', label='maximum reach')
    )
    legend_elements.append(
        Patch(facecolor='none', hatch='\\\\', edgecolor='k', label='IceCube arXiv:1709.03434')
    )

    legend = first_ax.legend(
        handles=legend_elements, prop=dict(size=11), loc='upper left',
        title='Excluded regions', framealpha=1., edgecolor='black',
        frameon=True
    )
    first_ax.set_zorder(10)
    plt.setp(legend.get_title(), fontsize='11')
    legend.get_frame().set_linestyle('-')

    if show_data: ybound = 0.65
    else: ybound = 0.73
    if args.data is DataType.REAL and show_data:
        # fig.text(0.295, 0.684, 'IceCube Preliminary', color='red', fontsize=13,
        fig.text(0.278, ybound, r'\bf IceCube Preliminary', color='red', fontsize=13,
                 ha='center', va='center', zorder=11)
    else:
        fig.text(0.278, ybound, r'\bf IceCube Simulation', color='red', fontsize=13,
                 ha='center', va='center', zorder=11)

    misc_utils.make_dir(outfile)
    for of in outformat:
        print 'Saving plot as {0}'.format(outfile+'.'+of)
        fig.savefig(outfile+'.'+of, bbox_inches='tight', dpi=150)


def plot_sens_fixed_angle(data, outfile, outformat, args):
    print 'Making FIXED_ANGLE sensitivity plot'

    colour = {0:'red', 1:'blue', 2:'green', 3:'purple', 4:'orange', 5:'black'}
    xticks = [r'$\mathcal{O}_{e\mu}$', r'$\mathcal{O}_{e\tau}$', r'$\mathcal{O}_{\mu\tau}$']
    argsc = deepcopy(args)

    LV_atmo_90pc_limits = {
        3: (2E-24, 1E-1),
        4: (2.7E-28, 3.16E-25),
        5: (1.5E-32, 1.12E-27),
        6: (9.1E-37, 2.82E-30),
        7: (3.6E-41, 1.77E-32),
        8: (1.4E-45, 1.00E-34)
    }
    ylims = {
        3 : (-28, -22), 4 : (-34, -25), 5 : (-42, -28), 6 : (-48, -33),
        7 : (-54, -37), 8 : (-61, -40)
    }

    for idim in xrange(len(data)):
        dim = args.dimensions[idim]
        print '= dim', dim
        argsc.dimension = dim

        # yranges = [np.inf, -np.inf]
        legend_handles = []

        fig = plt.figure(figsize=(7, 5))
        ax = fig.add_subplot(111)
        ax.set_xlim(0, len(xticks)+1)
        ax.set_xticklabels([''] + xticks + [''], fontsize=16)
        ax.set_xlabel(r'BSM operator angle', fontsize=16)
        ax.set_ylabel(r'${\mathrm {log}}_{10} \left (\Lambda^{-1}' + \
                      get_units(dim) +r'\right )$', fontsize=17)

        ax.tick_params(axis='y', labelsize=15)

        for isrc in xrange(len(data[idim])):
            src = args.source_ratios[isrc]
            argsc.source_ratio = src
            print '== src', src

            for ian in xrange(len(data[idim][isrc])):
                print '=== an', ian
                scales, statistic = data[idim][isrc][ian].T
                tck, u = interpolate.splprep([scales, statistic], s=0)
                scales, statistic = interpolate.splev(np.linspace(0, 1, 1000), tck)
                min_idx = np.argmin(scales)
                null = statistic[min_idx]
                if args.stat_method is StatCateg.BAYESIAN:
                    reduced_ev = -(statistic - null)
                    al = scales[reduced_ev > np.log(10**(3/2.))] # Strong degree of belief
                elif args.stat_method is StatCateg.FREQUENTIST:
                    reduced_ev = -2*(statistic - null)
                    al = scales[reduced_ev > 2.71] # 90% CL for 1 DOF via Wilks
                if len(al) == 0:
                    print 'No points for DIM {0} FRS {1}!'.format(dim, src)
                    continue
                if reduced_ev[-1] < np.log(10**(3/2.)) - 0.1:
                    print 'Peaked contour does not exclude large scales! For ' \
                        'DIM {0} FRS{1}!'.format(dim, src)
                    continue
                arr_len = 1.5
                lim = al[0]
                print 'limit = {0}'.format(lim)
                label = '{0} : {1} : {2}'.format(*misc_utils.solve_ratio(src))
                # if lim < yranges[0]: yranges[0] = lim-arr_len
                # if lim > yranges[1]: yranges[1] = lim+arr_len+2
                # if lim > yranges[1]: yranges[1] = lim
                xoff = 0.15
                line = plt.Line2D(
                    (ian+1-xoff, ian+1+xoff), (lim, lim), lw=2.5, color=colour[isrc], label=label
                )
                ax.add_line(line)
                if len(legend_handles) < isrc+1:
                    legend_handles.append(line)
                x_offset = isrc*xoff/2. - xoff/2.
                ax.annotate(
                    s='', xy=(ian+1+x_offset, lim+0.02), xytext=(ian+1+x_offset, lim-arr_len),
                    arrowprops={'arrowstyle': '<-', 'lw': 1.5, 'color':colour[isrc]}
                )
                if ian == 2:
                    lim = np.log10(LV_atmo_90pc_limits[dim][0])
                    ax.add_patch(patches.Rectangle(
                        (ian+1-xoff, lim), 2*xoff, 100, fill=True,
                        facecolor='orange', alpha=0.3, linewidth=1, edgecolor='k'
                    ))
                    ax.annotate(s='IC atmospheric', xy=(ian+1, lim+0.13),
                                color='red', rotation=90, fontsize=7,
                                horizontalalignment='center',
                                verticalalignment='bottom')

        # try:
        #     yranges = (myround(yranges[0], up=True), myround(yranges[1], down=True))
        #     ax.set_ylim(yranges)
        # except: pass
        ax.set_ylim(ylims[dim])

        legend = ax.legend(handles=legend_handles, prop=dict(size=10), loc='lower left',
                           title=r'$\nu_e:\nu_\mu:\nu_\tau{\rm\:\:at\:\:source}$',
                           framealpha=1., edgecolor='black')
        plt.setp(legend.get_title(), fontsize='10')
        legend.get_frame().set_linestyle('-')

        an_text = 'Dimension {0}'.format(dim)
        at = AnchoredText(
            an_text, prop=dict(size=10), frameon=True, loc=2
        )
        at.patch.set_boxstyle("round,pad=0.1,rounding_size=0.5")
        ax.add_artist(at)

        fig.text(0.52, 0.8, r'Excluded', color='red', fontsize=16, ha='center',
                 va='center', fontweight='bold')
        # fig.text(0.52, 0.76, r'with strong evidence', color='red', fontsize=11,
        #          ha='center', va='center')

        if args.data is DataType.REAL:
            fig.text(0.77, 0.14, 'IceCube Preliminary', color='red', fontsize=10,
                     ha='center', va='center')
        elif args.data in [DataType.ASIMOV, DataType.REALISATION]:
            fig.text(0.77, 0.14, 'IceCube Simulation', color='red', fontsize=10,
                     ha='center', va='center')

        for ymaj in ax.yaxis.get_majorticklocs():
            ax.axhline(y=ymaj, ls=':', color='gray', alpha=0.2, linewidth=1)
        for xmaj in ax.xaxis.get_majorticklocs():
            ax.axvline(x=xmaj, ls=':', color='gray', alpha=0.2, linewidth=1)

        out = outfile + '_DIM{0}'.format(dim)
        misc_utils.make_dir(out)
        for of in outformat:
            print 'Saving plot as {0}'.format(out+'.'+of)
            fig.savefig(out+'.'+of, bbox_inches='tight', dpi=150)


def plot_sens_corr_angle(data, outfile, outformat, args):
    print 'Making CORR_ANGLE sensitivity plot'

    labels = [r'$sin^2(2\mathcal{O}_{12})$',
              r'$sin^2(2\mathcal{O}_{13})$',
              r'$sin^2(2\mathcal{O}_{23})$']

    argsc = deepcopy(args)
    for idim in xrange(len(data)):
        dim = args.dimensions[idim]
        print '= dim', dim
        argsc.dimension = dim
        for isrc in xrange(len(data[idim])):
            src = args.source_ratios[isrc]
            argsc.source_ratio = src
            print '== src', src
            for ian in xrange(len(data[idim][isrc])):
                print '=== an', ian

                d = data[idim][isrc][ian].reshape(
                    len(data[idim][isrc][ian])**2, 3
                )

                fig = plt.figure(figsize=(7, 5))
                ax = fig.add_subplot(111)
                ax.set_ylim(0, 1)
                ax.set_ylabel(labels[ian])
                ax.set_xlabel(r'${\rm log}_{10} \left (\Lambda^{-1}'+get_units(dim)+r'\right )$')

                xranges = [np.inf, -np.inf]
                legend_handles = []

                y = d[:,0]
                x = d[:,1]
                z = d[:,2]

                print 'x', x
                print 'y', y
                print 'z', z
                null_idx = np.argmin(x)
                null = z[null_idx]
                print 'null = {0}, for scale = {1}'.format(null, x[null_idx])

                if args.stat_method is StatCateg.BAYESIAN:
                    z_r = -(z - null)
                elif args.stat_method is StatCateg.FREQUENTIST:
                    z_r = -2*(z - null)
                print 'scale', x
                print 'reduced ev', z_r

                pdat = np.array([x, y, z_r, np.ones(x.shape)]).T
                print 'pdat', pdat
                pdat_clean = []
                for d in pdat:
                    if not np.any(np.isnan(d)): pdat_clean.append(d)
                pdat = np.vstack(pdat_clean)
                sort_column = 3
                pdat_sorted = pdat[pdat[:,sort_column].argsort()]
                uni, c = np.unique(pdat[:,sort_column], return_counts=True)
                print uni, c
                print len(uni)
                print np.unique(c)
                
                n = len(uni)
                assert len(np.unique(c)) == 1
                c = c[0]
                col_array = []
                for col in pdat_sorted.T:
                    col_array.append(col.reshape(n, c))
                col_array = np.stack(col_array)
                sep_arrays = []
                for x_i in xrange(n):
                    sep_arrays.append(col_array[:,x_i])
                
                print len(sep_arrays)
                sep_arrays = sep_arrays[0][:3]
                print sep_arrays

                if args.stat_method is StatCateg.BAYESIAN:
                    reduced_pdat_mask = (sep_arrays[2] > np.log(10**(3/2.))) # Strong degree of belief
                elif args.stat_method is StatCateg.FREQUENTIST:
                    reduced_pdat_mask = (sep_arrays[2] > 4.61) # 90% CL for 2 DOFS via Wilks
                reduced_pdat = sep_arrays.T[reduced_pdat_mask].T
                print 'reduced_pdat', reduced_pdat

                ax.tick_params(axis='x', labelsize=11)
                ax.tick_params(axis='y', labelsize=11)

                mini, maxi = np.min(x), np.max(x)
                ax.set_xlim((mini, maxi))
                ax.set_ylim(0, 1)
                ax.grid(b=False)
                
                x_v = reduced_pdat[0].round(decimals=4)
                y_v = reduced_pdat[1].round(decimals=4)
                uniques = np.unique(x_v)
                print 'uniques', uniques
                if len(uniques) == 1: continue
                bw = np.min(np.diff(uniques))
                print 'bw', bw
                print np.diff(uniques)
                uni_x_split = np.split(uniques, np.where(np.diff(uniques) > bw*(1e20))[0] + 1)
                print 'len(uni_x_split)', len(uni_x_split)
                for uni_x in uni_x_split:
                    p_x_l, p_y_l = [], []
                    p_x_u, p_y_u = [], []
                    for uni in uni_x:
                        idxes = np.where(x_v == uni)[0]
                        ymin, ymax = 1, 0
                        for idx in idxes:
                            if y_v[idx] < ymin: ymin = y_v[idx]
                            if y_v[idx] > ymax: ymax = y_v[idx]
                        p_x_l.append(uni)
                        p_y_l.append(ymin)
                        p_x_u.append(uni)
                        p_y_u.append(ymax)
                    p_x_l, p_y_l = np.array(p_x_l, dtype=np.float64), np.array(p_y_l, dtype=np.float64)
                    p_x_u, p_y_u = np.array(list(reversed(p_x_u)), dtype=np.float64), np.array(list(reversed(p_y_u)), dtype=np.float64)
                    p_x = np.hstack([p_x_l, p_x_u])
                    p_y = np.hstack([p_y_l, p_y_u])
                    p_x = np.r_[p_x, p_x[0]]
                    p_y = np.r_[p_y, p_y[0]]
                    print 'p_x', p_x
                    print 'p_y', p_y
                    try:
                        tck, u = interpolate.splprep([p_x, p_y], s=1e-5, per=True)
                        xi, yi = interpolate.splev(np.linspace(0, 1, 1000), tck)
                    except:
                        xi, yi = p_x, p_y
                    ax.fill(xi, yi, 'r', edgecolor='r', linewidth=1)

                out = outfile + '_DIM{0}_SRC{1}_AN{2}'.format(dim, isrc, ian)
                misc_utils.make_dir(out)
                for of in outformat:
                    print 'Saving plot as {0}'.format(out+'.'+of)
                    fig.savefig(out+'.'+of, bbox_inches='tight', dpi=150)
