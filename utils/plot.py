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
import numpy.ma as ma
from scipy.interpolate import splev, splprep
from scipy.ndimage.filters import gaussian_filter

import matplotlib as mpl
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
mpl.use('Agg')

from matplotlib import rc
from matplotlib import pyplot as plt
from matplotlib.offsetbox import AnchoredText
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from getdist import plots, mcsamples

import ternary
from ternary.heatmapping import polygon_generator

import shapely.geometry as geometry

from utils import misc as misc_utils
from utils.enums import DataType, EnergyDependance, str_enum
from utils.enums import Likelihood, ParamTag, StatCateg, Texture
from utils.fr import angles_to_u, angles_to_fr


BAYES_K = 1.   # Substantial degree of belief.
# BAYES_K = 3/2. # Strong degree of belief.
# BAYES_K = 2.   # Very strong degree of belief
# BAYES_K = 5/2. # Decisive degree of belief


if os.path.isfile('./plot_llh/paper.mplstyle'):
    plt.style.use('./plot_llh/paper.mplstyle')
elif os.environ.get('GOLEMSOURCEPATH') is not None:
    plt.style.use(os.environ['GOLEMSOURCEPATH']+'/GolemFit/scripts/paper/paper.mplstyle')
if 'submitter' in socket.gethostname():
    rc('text', usetex=False)


def gen_figtext(args):
    """Generate the figure text."""
    t = ''
    t += 'Source flavour ratio = [{0}]'.format(solve_ratio(args.source_ratio))
    if args.data in [DataType.ASIMOV, DataType.REALISATION]:
        t += '\nIC injected flavour ratio = [{0}]'.format(
            solve_ratio(args.injected_ratio)
        )
    t += '\nDimension = {0}'.format(args.dimension)
    return t


def texture_label(x):
    if x == Texture.OEU:
        return r'$\mathcal{O}_{e\mu}$'
    elif x == Texture.OET:
        return r'$\mathcal{O}_{e\tau}$'
    elif x == Texture.OUT:
        return r'$\mathcal{O}_{\mu\tau}$'
    else:
        raise AssertionError


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


def chainer_plot(infile, outfile, outformat, args, llh_paramset, fig_text=None):
    """Make the triangle plot."""
    if hasattr(args, 'plot_elements'):
        if not args.plot_angles and not args.plot_elements:
            return
    elif not args.plot_angles:
        return

    if not isinstance(infile, np.ndarray):
        raw = np.load(infile)
    else:
        raw = infile
    print 'raw.shape', raw.shape
    print 'raw', raw

    misc_utils.make_dir(outfile)
    if fig_text is None:
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
            plt.text(0.8, 0.9, 'IceCube Preliminary', color='red', fontsize=15,
                     ha='center', va='center')
        elif args.data in [DataType.ASIMOV, DataType.REALISATION]:
            plt.text(0.8, 0.9, 'IceCube Simulation', color='red', fontsize=15,
                     ha='center', va='center')

        for of in outformat:
            print 'Saving', outfile+'_angles.'+of
            g.export(outfile+'_angles.'+of)

    if not hasattr(args, 'plot_elements'):
        return

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
        if args.fix_mixing is not MixingScenario.NONE:
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
            plt.text(0.8, 0.7, 'IceCube Preliminary', color='red', fontsize=15,
                     ha='center', va='center')
        elif args.data in [DataType.ASIMOV, DataType.REALISATION]:
            plt.text(0.8, 0.7, 'IceCube Simulation', color='red', fontsize=15,
                     ha='center', va='center')

        mpl.pyplot.figtext(0.5, 0.7, fig_text, fontsize=15)
        for of in outformat:
            print 'Saving', outfile+'_elements'+of
            g.export(outfile+'_elements.'+of)


def plot_statistic(data, outfile, outformat, args, scale_param, label=None):
    """Make MultiNest factor or LLH value plot."""
    print 'Making Statistic plot'
    fig_text = gen_figtext(args)
    if label is not None: fig_text += '\n' + label

    print 'data', data
    print 'data.shape', data.shape
    scales, statistic = ma.compress_rows(data).T
    try:
        tck, u = splprep([scales, statistic], s=0)
    except:
        return
    sc, st = splev(np.linspace(0, 1, 10000), tck)
    scales_rm = sc[sc >= scales[1]]
    statistic_rm = st[sc >= scales[1]]

    min_idx = np.argmin(scales)
    null = statistic[min_idx]
    fig_text += '\nnull lnZ = {0:.2f}'.format(null)

    if args.stat_method is StatCateg.BAYESIAN:
        reduced_ev = -(statistic_rm - null)
    elif args.stat_method is StatCateg.FREQUENTIST:
        reduced_ev = -2*(statistic_rm - null)

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)

    ax.set_xlim(np.log10(args.scale_region))
    ax.set_xlabel(r'${\mathrm {log}}_{10} \left (\Lambda^{-1}' + \
                  get_units(args.dimension) +r'\right )$', fontsize=16)
    if args.stat_method is StatCateg.BAYESIAN:
        ax.set_ylabel(r'log(Bayes Factor)')
    elif args.stat_method is StatCateg.FREQUENTIST:
        ax.set_ylabel(r'$-2\Delta {\rm LLH}$')

    # ymin = np.round(np.min(reduced_ev) - 1.5)
    # ymax = np.round(np.max(reduced_ev) + 1.5)
    # ax.set_ylim((ymin, ymax))

    ax.plot(scales_rm, reduced_ev)

    ax.axhline(y=np.log(10**(BAYES_K)), color='red', alpha=1., linewidth=1.3)

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
                al = scales[reduced_ev > np.log(10**(BAYES_K))] # Strong degree of belief
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


def plot_table_sens(data, outfile, outformat, args):
    print 'Making TABLE sensitivity plot'
    argsc = deepcopy(args)

    dims = args.dimensions
    srcs = args.source_ratios
    if args.texture is Texture.NONE:
        textures = [Texture.OEU, Texture.OET, Texture.OUT]
    else:
        textures = [args.texture]

    if len(srcs) > 3:
        raise NotImplementedError

    xlims = (-60, -20)
    ylims = (0.5, 1.5)

    LV_atmo_90pc_limits = {
        3: (2E-24, 1E-1),
        4: (2.7E-28, 3.16E-25),
        5: (1.5E-32, 1.12E-27),
        6: (9.1E-37, 2.82E-30),
        7: (3.6E-41, 1.77E-32),
        8: (1.4E-45, 1.00E-34)
    }

    PS = 8.203e-20 # GeV^{-1}
    planck_scale = {
        5: PS,
        6: PS**2,
        7: PS**3,
        8: PS**4
    }

    colour = {0:'red', 1:'blue', 2:'green'}
    rgb_co = {0:[1,0,0], 1:[0,0,1], 2:[0.0, 0.5019607843137255, 0.0]}

    fig = plt.figure(figsize=(8, 6))
    gs = gridspec.GridSpec(dims, 1)
    gs.update(hspace=0.15)

    first_ax = None
    legend_log = []
    legend_elements = []

    for idim, dim in enumerate(dimensions):
        print '== dim', dim
        argsc.dimension = dim
        gs0 = gridspec.GridSpecFromSubplotSpec(
            len(textures), 1, subplot_spec=gs[idim], hspace=0
        )

        for itex, tex in enumerate(textures):
            argcs.texture = tex
            ylabel = texture_label(texture)
            # if angles == 2 and ian == 0: continue
            ax = fig.add_subplot(gs0[itex])

            if first_ax is None:
                first_ax = ax

            ax.set_xlim(xlims)
            ax.set_ylim(ylims)
            ax.set_yticks([1.])
            ax.set_yticklabels([ylabel], fontsize=13)
            ax.yaxis.tick_right()
            for xmaj in ax.xaxis.get_majorticklocs():
                ax.axvline(x=xmaj, ls=':', color='gray', alpha=0.2, linewidth=1)
            ax.get_xaxis().set_visible(False)
            # TODO(shivesh): check this
            if itex == len(tex) - 2:
                ax.spines['bottom'].set_alpha(0.6)
            elif itex == len(tex) - 1:
                ax.text(
                    -0.04, ylim[0], '$d = {0}$'.format(dim), fontsize=16,
                    rotation='90', verticalalignment='center',
                    transform=ax.transAxes
                )
                dim_label_flag = False
                ax.spines['top'].set_alpha(0.6)
                ax.spines['bottom'].set_alpha(0.6)

            for isrc, src in enumerate(srcs):
                print '== src', src
                argsc.source_ratio = src

                if dim in planck_scale:
                    ps = np.log10(planck_scale[dim])
                    if ps < xlims[0]:
                        ax.annotate(
                            s='', xy=(xlims[0], 1), xytext=(xlims[0]+1, 1),
                            arrowprops={'arrowstyle': '->, head_length=0.2',
                                        'lw': 1, 'color':'purple'}
                        )
                    elif ps > xlims[1]:
                        ax.annotate(
                            s='', xy=(xlims[1]-1, 1), xytext=(xlims[1], 1),
                            arrowprops={'arrowstyle': '<-, head_length=0.2',
                                        'lw': 1, 'color':'purple'}
                        )
                    else:
                        ax.axvline(x=ps, color='purple', alpha=1., linewidth=1.5)

                scales, statistic = ma.compress_rows(data[idim][isrc][itex]).T
                try:
                    tck, u = splprep([scales, statistic], s=0)
                except:
                    return
                sc, st = splev(np.linspace(0, 1, 10000), tck)
                scales_rm = sc[sc >= scales[1]]
                statistic_rm = st[sc >= scales[1]]

                min_idx = np.argmin(scales)
                null = statistic[min_idx]
                if args.stat_method is StatCateg.BAYESIAN:
                    reduced_ev = -(statistic_rm - null)
                    al = scales_rm[reduced_ev > np.log(10**(BAYES_K))]
                elif args.stat_method is StatCateg.FREQUENTIST:
                    reduced_ev = -2*(statistic_rm - null)
                    al = scales_rm[reduced_ev > 2.71] # 90% CL for 1 DOF via Wilks
                if len(al) == 0:
                    print 'No points for DIM {0} FRS {1}!'.format(dim, src)
                    continue
                if reduced_ev[-1] < np.log(10**(BAYES_K)) - 0.1:
                    print 'Peaked contour does not exclude large scales! For ' \
                        'DIM {0} FRS{1}!'.format(dim, src)
                    continue
                lim = al[0]
                print 'limit = {0}'.format(lim)

                ax.axvline(x=lim, color=colour[isrc], alpha=1., linewidth=1.5)
                ax.add_patch(patches.Rectangle(
                    (lim, ylim[0]), 100, np.diff(ylim), fill=True, facecolor=colour[isrc],
                    alpha=0.3, linewidth=0
                ))

                if isrc not in legend_log:
                    legend_log.append(isrc)
                    label = '{0} at source'.format(misc_utils.solve_ratio(src))
                    legend_elements.append(
                        Patch(facecolor=rgb_co[isrc]+[0.3],
                              edgecolor=rgb_co[isrc]+[1], label=label)
                    )

            if itex == 2:
                LV_lim = np.log10(LV_atmo_90pc_limits[dim])
                ax.add_patch(patches.Rectangle(
                    (LV_lim[1], ylim[0]), LV_lim[0]-LV_lim[1], np.diff(ylim),
                    fill=False, hatch='\\\\'
                ))

    ax.get_xaxis().set_visible(True)
    ax.set_xlabel(r'${\rm New\:Physics\:Scale}\:[\:{\rm log}_{10} (\Lambda^{-1}\:/\:{\rm GeV}^{-d+4})\: ]$',
                 fontsize=19)
    ax.tick_params(axis='x', labelsize=16)

    purple = [0.5019607843137255, 0.0, 0.5019607843137255]
    legend_elements.append(
        Patch(facecolor=purple+[0.7], edgecolor=purple+[1], label='Planck Scale Expectation')
    )
    legend_elements.append(
        Patch(facecolor='none', hatch='\\\\', edgecolor='k', label='IceCube, Nature.Phy.14,961(2018)')
    )

    legend = first_ax.legend(
        handles=legend_elements, prop=dict(size=11), loc='upper left',
        title='Excluded regions', framealpha=1., edgecolor='black',
        frameon=True
    )
    first_ax.set_zorder(10)
    plt.setp(legend.get_title(), fontsize='11')
    legend.get_frame().set_linestyle('-')

    if show_data: ybound = 0.595
    else: ybound = 0.73
    if args.data is DataType.REAL and show_data:
        # fig.text(0.295, 0.684, 'IceCube Preliminary', color='red', fontsize=13,
        fig.text(0.278, ybound, r'\bf IceCube Preliminary', color='red', fontsize=13,
                 ha='center', va='center', zorder=11)
    elif args.data is DataType.REALISATION:
        fig.text(0.278, ybound-0.05, r'\bf IceCube Simulation', color='red', fontsize=13,
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
                tck, u = splprep([scales, statistic], s=0)
                scales, statistic = splev(np.linspace(0, 1, 1000), tck)
                min_idx = np.argmin(scales)
                null = statistic[min_idx]
                if args.stat_method is StatCateg.BAYESIAN:
                    reduced_ev = -(statistic - null)
                    al = scales[reduced_ev > np.log(10**(BAYES_K))] # Strong degree of belief
                elif args.stat_method is StatCateg.FREQUENTIST:
                    reduced_ev = -2*(statistic - null)
                    al = scales[reduced_ev > 2.71] # 90% CL for 1 DOF via Wilks
                if len(al) == 0:
                    print 'No points for DIM {0} FRS {1}!'.format(dim, src)
                    continue
                if reduced_ev[-1] < np.log(10**(BAYES_K)) - 0.1:
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
            fig_text = gen_figtext(argsc)
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
                    reduced_pdat_mask = (sep_arrays[2] > np.log(10**(BAYES_K))) # Strong degree of belief
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
                        tck, u = splprep([p_x, p_y], s=1e-5, per=True)
                        xi, yi = splev(np.linspace(0, 1, 1000), tck)
                    except:
                        xi, yi = p_x, p_y
                    ax.fill(xi, yi, 'r', edgecolor='r', linewidth=1)

                at = AnchoredText(
                    fig_text, prop=dict(size=10), frameon=True, loc=4
                )
                at.patch.set_boxstyle("round,pad=0.1,rounding_size=0.5")
                ax.add_artist(at)
                out = outfile + '_DIM{0}_SRC{1}_AN{2}'.format(dim, isrc, ian)
                misc_utils.make_dir(out)
                for of in outformat:
                    print 'Saving plot as {0}'.format(out+'.'+of)
                    fig.savefig(out+'.'+of, bbox_inches='tight', dpi=150)


def cmap_discretize(cmap, N):
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in xrange(N+1) ]
    # Return colormap object.
    return mpl.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)


def get_tax(ax, scale, ax_labels):
    ax.set_aspect('equal')

    # Boundary and Gridlines
    fig, tax = ternary.figure(ax=ax, scale=scale)

    # Draw Boundary and Gridlines
    tax.boundary(linewidth=2.0)
    tax.gridlines(color='grey', multiple=scale/5., linewidth=1.0, alpha=0.4, ls='--')
    tax.gridlines(color='grey', multiple=scale/10., linewidth=0.5, alpha=0.4, ls=':')

    # Set Axis labels and Title
    fontsize = 23
    tax.bottom_axis_label(ax_labels[0], fontsize=fontsize+8, position=(0.55, -0.20/2, 0.5), rotation=0)
    tax.right_axis_label(ax_labels[1], fontsize=fontsize+8, offset=0.2, rotation=0)
    tax.left_axis_label(ax_labels[2], fontsize=fontsize+8, offset=0.2, rotation=0)

    # Remove default Matplotlib axis
    tax.get_axes().axis('off')
    tax.clear_matplotlib_ticks()

    # Set ticks
    # ticks = np.linspace(0, 1, 6)
    # tax.ticks(ticks=ticks, locations=ticks*scale, axis='blr', linewidth=1,
    #           offset=0.03, fontsize=fontsize, tick_formats='%.1f')
    tax.ticks()

    tax._redraw_labels()

    return tax

def triangle_project(frs, llh, outfile, outformat, args, llh_paramset, fig_text):
    print 'Making triangle projection'
    fontsize = 23

    def alp(x):
        y = list(x)
        y[-1] = 0.4
        return y

    cmap = plt.get_cmap('jet', 10)
    cmap_g = cmap_discretize(
        mpl.colors.LinearSegmentedColormap.from_list(
            "", ["lime", "gold", "darkorange"]),
        10
    )
    cmap_b = cmap_discretize(
        mpl.colors.LinearSegmentedColormap.from_list(
            "", ["blue", "fuchsia", "darkmagenta"]),
        10
    )

    mean = np.mean(llh)
    sig = np.std(llh)
    max_scale = np.max(llh)
    min_scale = np.min(mean-sig)
    norm = mpl.colors.Normalize(vmin=min_scale, vmax=max_scale)

    colour = map(alp, map(cmap, map(norm, llh)))
    # colour = map(alp, map(cmap_g, map(norm, llh)))
    # colour = map(alp, map(cmap_b, map(norm, llh)))

    # Figure
    fig = plt.figure(figsize=(8, 8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[40, 1])
    gs.update(hspace=0.3, wspace=0.15)

    ax = fig.add_subplot(gs[0])
    tax = get_tax(ax, scale=1)

    # Plot
    tax.scatter(frs, marker='o', s=0.1, color=colour)

    # Contour TODO(shivesh)
    # tax.plot(contour_68_upper, linewidth=2.3, color='grey', zorder=0, alpha=0.6)
    # tax.plot(contour_68_lower, linewidth=2.3, color='grey', zorder=0, alpha=0.6)
    # tax.plot(contour_90_upper, linewidth=2.3, color='darkgrey', zorder=0, alpha=0.6)
    # tax.plot(contour_90_lower, linewidth=2.3, color='darkgrey', zorder=0, alpha=0.6)

    # Lines
    marker_style = dict(
        linestyle=' ', color='darkorange', linewidth=1.2,
        markersize=14, zorder=0
    )

    # Colorbar
    ax0 = fig.add_subplot(gs[1])
    cb = mpl.colorbar.ColorbarBase(ax0, cmap=cmap, norm=norm,
                                   orientation='vertical')
    cb.ax.tick_params(labelsize=fontsize-5)
    cb.set_label(r'$LLH$', fontsize=fontsize+5, labelpad=20,
                 horizontalalignment='center')

    misc_utils.make_dir(outfile)
    for of in outformat:
        print 'Saving plot as {0}'.format(outfile+'.'+of)
        fig.savefig(outfile+'.'+of, bbox_inches='tight', dpi=150)


def heatmap(data, scale, vmin=None, vmax=None, style='triangular'):
    for k, v in data.items():
        data[k] = np.array(v)
    style = style.lower()[0]
    if style not in ["t", "h", 'd']:
        raise ValueError("Heatmap style must be 'triangular', 'dual-triangular', or 'hexagonal'")

    vertices_values = polygon_generator(data, scale, style)

    vertices = []
    for v, value in vertices_values:
        vertices.append(v)
    return vertices


def flavour_contour(frs, ax, nbins, coverage, **kwargs):
    """Plot the flavour contour for a specified coverage."""
    # Histogram in flavour space
    H, b = np.histogramdd(
        (frs[:,0], frs[:,1], frs[:,2]),
        bins=(nbins+1, nbins+1, nbins+1), range=((0, 1), (0, 1), (0, 1))
    )
    H = H / np.sum(H)

    # 3D smoothing
    smoothing = 0.05
    H_s = gaussian_filter(H, sigma=smoothing)

    # Finding coverage
    H_r = np.ravel(H_s)
    H_rs = np.argsort(H_r)[::-1]
    H_crs = np.cumsum(H_r[H_rs])
    thres = np.searchsorted(H_crs, coverage/100.)
    mask_r = np.zeros(H_r.shape)
    mask_r[H_rs[:thres]] = 1
    mask = mask_r.reshape(H_s.shape)

    # Get vertices inside covered region
    binx = np.linspace(0, 1, nbins+1)
    interp_dict = {}
    for i in xrange(len(binx)):
        for j in xrange(len(binx)):
            for k in xrange(len(binx)):
                if mask[i][j][k] == 1:
                    interp_dict[(i, j, k)] = H_s[i, j, k]
    vertices = np.array(heatmap(interp_dict, nbins))
    points = vertices.reshape((len(vertices)*3, 2))

    # Convex full to find points forming exterior bound
    pc = geometry.MultiPoint(points)
    polygon = pc.convex_hull
    ex_cor = np.array(list(polygon.exterior.coords))

    # Join points with a spline
    tck, u = splprep([ex_cor.T[0], ex_cor.T[1]], s=0., per=1, k=1)
    xi, yi = map(np.array, splev(np.linspace(0, 1, 300), tck))

    # Spline again to smooth
    tck, u = splprep([xi, yi], s=0.4, per=1, k=3)
    xi, yi = map(np.array, splev(np.linspace(0, 1, 300), tck))
    ev_polygon = np.dstack((xi, yi))[0]

    def project_toflavour(p):
        """Convert from cartesian to flavour space."""
        x, y = p
        b = y / (np.sqrt(3)/2.)
        a = x - b/2.
        return [a, b, nbins-a-b]

    # Remove points interpolated outside flavour triangle
    f_ev_polygon = np.array(map(project_toflavour, ev_polygon))
    xf, yf, zf = f_ev_polygon.T
    mask = np.array((xf < 0) | (yf < 0) | (zf < 0) | (xf > nbins) |
                    (yf > nbins) | (zf > nbins))
    ev_polygon = np.dstack((xi[~mask], yi[~mask]))[0]

    # Plot
    ax.plot(
        ev_polygon.T[0], ev_polygon.T[1], label=r'{0}\%'.format(int(coverage)),
        **kwargs
    )
    ax.scatter(points.T[0], points.T[1], marker='o', s=2, alpha=1, color=color,
               zorder=3)

def plot_source_ternary(data, outfile, outformat, args):
    """Ternary plot in the source flavour space for each operator texture."""
    r_data = np.full((data.shape[0], data.shape[2], data.shape[1], data.shape[3], data.shape[4]), np.nan)
    for idim in xrange(data.shape[0]):
        for isrc in xrange(data.shape[1]):
            for isce in xrange(data.shape[2]):
                r_data[idim][isce][isrc] = data[idim][isrc][isce]
    r_data = ma.masked_invalid(r_data)
    print r_data.shape, 'r_data.shape'
    nsrcs = int(len(args.source_ratios) / 3)

    for idim, dim in enumerate(args.dimensions):
        print '|||| DIM = {0}, {1}'.format(idim, dim)
        for isce in xrange(r_data.shape[1]):
            print '|||| SCEN = {0}'.format(str_enum(MixingScenario(isce+1)))
            fig = plt.figure(figsize=(8, 8))
            ax  = fig.add_subplot(111)
            tax = get_tax(ax, scale=nsrcs)
            interp_dict = {}
            for isrc, src in enumerate(args.source_ratios):
                src_sc = tuple(np.array(src)*(nsrcs-1))
                print '|||| SRC = {0}'.format(src)
                scales, statistic = ma.compress_rows(r_data[idim][isce][isrc]).T
                print 'scales', scales
                print 'statistic', statistic

                try:
                    tck, u = splprep([scales, statistic], s=0)
                except:
                    interp_dict[src_sc] = -60
                    continue
                # sc, st = splev(np.linspace(0, 1, 10000), tck)
                sc, st = splev(np.linspace(0, 1, 20), tck)
                scales_rm = sc[sc >= scales[1]]
                statistic_rm = st[sc >= scales[1]]
                # scales_rm = sc
                # statistic_rm = st

                min_idx = np.argmin(scales)
                null = statistic[min_idx]
                if args.stat_method is StatCateg.BAYESIAN:
                    reduced_ev = -(statistic_rm - null)
                    print 'reduced_ev', reduced_ev
                    al = scales_rm[reduced_ev > np.log(10**(BAYES_K))]
                else:
                    assert 0
                if len(al) == 0:
                    print 'No points for DIM {0} FRS {1}!'.format(dim, src)
                    interp_dict[src_sc] = -60
                    continue
                if reduced_ev[-1] < np.log(10**(BAYES_K)) - 0.1:
                    print 'Peaked contour does not exclude large scales! For ' \
                        'DIM {0} FRS{1}!'.format(dim, src)
                    interp_dict[src_sc] = -60
                    continue
                lim = al[0]
                print 'limit = {0}'.format(lim)

                interp_dict[src_sc] = lim
            print 'interp_dict', interp_dict
            print
            print 'vertices', heatmap(interp_dict, nsrcs)
            print
            tax.heatmap(interp_dict, scale=nsrcs, vmin=-60, vmax=-30)
            misc_utils.make_dir(outfile)
            for of in outformat:
                print 'Saving plot as {0}'.format(outfile+'_SCEN{0}.'.format(isce)+of)
                fig.savefig(outfile+'_SCEN{0}.'.format(isce)+of, bbox_inches='tight', dpi=150)
            print 'nsrcs', nsrcs
            assert 0


def plot_x(data, outfile, outformat, args):
    """Limit plot as a function of the source flavour ratio for each operator
    texture."""
    print 'Making X sensitivity plot'
    dims = args.dimensions
    srcs = args.source_ratios
    x_arr = [i[0] for i in srcs]
    if args.texture is Texture.NONE:
        textures = [Texture.OEU, Texture.OET, Texture.OUT]
    else:
        textures = [args.texture]

    # Rearrange data structure
    r_data = np.full((
        data.shape[0], data.shape[2], data.shape[1], data.shape[3], data.shape[4]
    ), np.nan)

    for idim in xrange(data.shape[0]):
        for isrc in xrange(data.shape[1]):
            for itex in xrange(data.shape[2]):
                r_data[idim][itex][isrc] = data[idim][isrc][itex]
    r_data = ma.masked_invalid(r_data)
    print r_data.shape, 'r_data.shape'

    for idim, dim in enumerate(dims):
        print '|||| DIM = {0}, {1}'.format(idim, dim)
        fig = plt.figure(figsize=(4, 4))
        ax  = fig.add_subplot(111)
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)
        ax.set_xlabel(r'$x$', fontsize=18)
        ax.set_ylabel(r'${\rm New\:Physics\:Scale}\:[\:{\rm log}_{10} (\Lambda_{d='+str(dim)+r'}^{-1}\:/\:'+get_units(args.dimension)+r')\: ]$', fontsize=18)
        ax.set_xlim(0, 1)
        for itex, tex in enumerate(texture):
            print '|||| TEX = {0}'.format(texture)
            lims = np.full(len(srcs), np.nan)
            for isrc, src in enumerate(srcs):
                x = src[0]
                print '|||| X = {0}'.format(x)
                scales, statistic = ma.compress_rows(r_data[idim][itex][isrc]).T
                print 'scales', scales
                print 'statistic', statistic

                try:
                    tck, u = splprep([scales, statistic], s=0)
                except:
                    continue
                sc, st = splev(np.linspace(0, 1, 10000), tck)
                scales_rm = sc[sc >= scales[1]]
                statistic_rm = st[sc >= scales[1]]

                min_idx = np.argmin(scales)
                null = statistic[min_idx]
                if args.stat_method is StatCateg.BAYESIAN:
                    reduced_ev = -(statistic_rm - null)
                    print 'reduced_ev', reduced_ev
                    al = scales_rm[reduced_ev > np.log(10**(BAYES_K))]
                else:
                    assert 0
                if len(al) == 0:
                    print 'No points for DIM {0} X {1}!'.format(dim, x)
                    continue
                if reduced_ev[-1] < np.log(10**(BAYES_K)) - 0.1:
                    print 'Peaked contour does not exclude large scales! For ' \
                        'DIM {0} X {1}!'.format(dim, x)
                    continue
                lim = al[0]
                print 'limit = {0}'.format(lim)

                lims[isrc] = lim
            lims = ma.masked_invalid(lims)
            print 'lims', lims
            ax.scatter(x_arr, lims)
        for ymaj in ax.yaxis.get_majorticklocs():
            ax.axhline(y=ymaj, ls=':', color='gray', alpha=0.3, linewidth=1)
        for xmaj in be:
            ax.axvline(x=xmaj, ls=':', color='gray', alpha=0.3, linewidth=1)
        ax.legend()
        misc_utils.make_dir(outfile)
        for of in outformat:
            print 'Saving plot as {0}'.format(outfile+'_DIM{0}.'.format(dim)+of)
            fig.savefig(outfile+'_DIM{0}.'.format(dim)+of, bbox_inches='tight', dpi=150)
