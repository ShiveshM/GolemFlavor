#! /usr/bin/env python
# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : April 14, 2018

"""
HESE BSM flavour ratio sensivity plotting script
"""

from __future__ import absolute_import, division

import os

import numpy as np
import numpy.ma as ma

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import rc
from matplotlib import pyplot as plt
from matplotlib.offsetbox import AnchoredText

import fr
from utils import misc as misc_utils
from utils.fr import normalise_fr
from utils.plot import myround, get_units


rc('text', usetex=False)
rc('font', **{'family':'serif', 'serif':['Computer Modern'], 'size':18})

fix_sfr_mfr = [
    (1, 1, 1, 1, 2, 0),
    # (1, 1, 1, 1, 0, 0),
    # (1, 1, 1, 0, 1, 0),
]

# FR
dimension         = [3, 4, 5, 6, 7, 8]
sigma_ratio       = ['0.01']
energy_dependance = 'spectral'
spectral_index    = -2
binning           = [1e4, 1e7, 5]
fix_mixing        = 'False'
fix_mixing_almost = 'False'
scale_region      = "1E10"

# Likelihood
likelihood = 'golemfit'

# Nuisance
convNorm        = 1.
promptNorm      = 0.
muonNorm        = 1.
astroNorm       = 6.9
astroDeltaGamma = 2.5

# GolemFit
ast  = 'p2_0'
data = 'real'

# Bayes Factor
bayes_bins        = 100
bayes_live_points = 1000
bayes_tolerance   = 0.01
bayes_eval_bin    = True # set to 'all' to run normally

# Plot
plot_bayes        = False
plot_angles_limit = False
plot_angles_corr  = True
outformat         = ['png']
significance      = np.log(10**(1/2.))
# significance      = np.log(10**(3/2.))


bayes_array = ma.masked_equal(np.zeros((len(dimension), len(fix_sfr_mfr), bayes_bins, 2)), 0)
angles_lim_array = np.zeros((len(dimension), len(fix_sfr_mfr), 3, bayes_bins, 2))
for i_dim, dim in enumerate(dimension):
    if energy_dependance == 'mono':
        outchain_head = '/data/user/smandalia/flavour_ratio/data/{0}/DIM{1}/{2:.0E}'.format(likelihood, dim, en)
    elif energy_dependance == 'spectral':
        outchain_head = '/data/user/smandalia/flavour_ratio/data/{0}/DIM{1}/SI_{2}'.format(likelihood, dim, spectral_index)

    bayes_output = 'None'
    angles_lim_output = 'None'
    angles_corr_output = 'None'
    for sig in sigma_ratio:
        for i_frs, frs in enumerate(fix_sfr_mfr):
            outchains = outchain_head + '/fix_ifr/{0}/'.format(str(sig).replace('.', '_'))
            if plot_bayes:
                bayes_output = outchains + '/bayes_factor/'
            if plot_angles_limit:
                angles_lim_output = outchains + '/angles_limit/'
            if plot_angles_corr:
                angles_corr_output = outchains + '/angles_corr/'

            argstring = '--measured-ratio {0} {1} {2} --fix-source-ratio True --source-ratio {3} {4} {5} --dimension {6} --seed 24 --outfile {7} --run-mcmc False --likelihood {8} --plot-angles False --bayes-output {9} --angles-lim-output {10} --bayes-bins {11} --angles-corr-output'.format(frs[0], frs[1], frs[2], frs[3], frs[4], frs[5], dim, outchains, likelihood, bayes_output, angles_lim_output, bayes_bins, angles_corr_output)
            args = fr.parse_args(argstring)
            fr.process_args(args)
            # misc_utils.print_args(args)

            if plot_bayes:
                infile = args.bayes_output+'/fr_evidence'+misc_utils.gen_identifier(args)
            if plot_angles_limit:
                infile = args.angles_lim_output+'/fr_an_evidence'+misc_utils.gen_identifier(args)
            if plot_angles_corr:
                infile = args.angles_corr_output+'/fr_co_evidence' + misc_utils.gen_identifier(args)
            scan_scales = np.linspace(
                np.log10(args.scale_region[0]), np.log10(args.scale_region[1]), args.bayes_bins
            )
            raw = []
            fail = 0
            for sc in scan_scales:
                try:
                    infile_s = infile + '_scale_{0:.0E}'.format(np.power(10, sc))
                    lf = np.load(infile_s+'.npy')
                    if plot_angles_limit:
                        if len(lf.shape) == 3: lf = lf[:,0,:]
                    if plot_angles_corr:
                        # TODO(shivesh)
                        assert 0
                        if len(lf.shape) == 3: lf = lf[:,0,:]
                    raw.append(lf)
                except IOError:
                    fail += 1
                    print 'failed to open {0}'.format(infile_s)
                    if plot_bayes:
                        raw.append([0, 0])
                    if plot_angles_limit:
                        raw.append(np.zeros((3, 2)))
                    pass
            print 'failed to open {0} files'.format(fail)

            if plot_bayes:
                raw = np.vstack(raw)
            if plot_angles_limit:
                raw = np.vstack(raw).reshape(args.bayes_bins, 3, 2)
                a = ma.masked_equal(np.zeros((3, args.bayes_bins, 2)), 0)
                for i_x, x in enumerate(raw):
                    for i_y, y in enumerate(x):
                        a[i_y][i_x] = ma.masked_equal(y, 0)

            if plot_bayes:
                bayes_array[i_dim][i_frs] = ma.masked_equal(raw, 0)

            if plot_angles_limit:
                # TODO(shivesh)
                angles_lim_array[i_dim][i_frs] = ma.masked_equal(a, 0)

if plot_bayes:
    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)

    colour = {0:'red', 1:'blue', 2:'green', 3:'purple', 4:'orange', 5:'black'}
    yranges = [np.inf, -np.inf]
    legend_handles = []
    ax.set_xlim(dimension[0]-1, dimension[-1]+1)
    xticks = [''] + range(dimension[0], dimension[-1]+1) + ['']
    ax.set_xticklabels(xticks)
    ax.set_xlabel(r'BSM operator dimension ' + r'$d$')
    ax.set_ylabel(r'${\rm log}_{10} \Lambda / GeV^{-d+4}$')
    for i_dim, dim in enumerate(dimension):
        for i_frs, frs in enumerate(fix_sfr_mfr):
            scale, evidences = bayes_array[i_dim][i_frs].T
            null = evidences[np.argmin(scale)]
            # TODO(shivesh): negative or not?
            reduced_ev = -(evidences - null)
            al = scale[reduced_ev > significance]
            if len(al) > 0:
                label = '[{0}, {1}, {2}]'.format(frs[3], frs[4], frs[5])
                lim = al[0]
                print 'frs, dim, lim = ', frs, dim, lim
                if lim < yranges[0]: yranges[0] = lim
                if lim > yranges[1]: yranges[1] = lim+4
                line = plt.Line2D(
                    (dim-0.1, dim+0.1), (lim, lim), lw=3, color=colour[i_frs], label=label
                )
                ax.add_line(line)
                if i_dim == 0: legend_handles.append(line)
                x_offset = i_frs*0.05 - 0.05
                ax.annotate(
                    s='', xy=(dim+x_offset, lim), xytext=(dim+x_offset, lim+3),
                    arrowprops={'arrowstyle': '<-', 'lw': 1.2, 'color':colour[i_frs]}
                )

            else:
                print 'No points for DIM {0} FRS {1} NULL {2}!'.format(dim, frs, null)
    yranges = (myround(yranges[0], up=True), myround(yranges[1], down=True))
    ax.set_ylim(yranges)

    ax.legend(handles=legend_handles, prop=dict(size=8), loc='upper right')
    for ymaj in ax.yaxis.get_majorticklocs():
        ax.axhline(y=ymaj, ls=':', color='gray', alpha=0.4, linewidth=1)
    for xmaj in ax.xaxis.get_majorticklocs():
        ax.axvline(x=xmaj, ls=':', color='gray', alpha=0.4, linewidth=1)

    for of in outformat:
        fig.savefig('./images/bayes_factor.'+of, bbox_inches='tight', dpi=150)

if plot_angles_limit:
    colour = {0:'red', 1:'blue', 2:'green', 3:'purple', 4:'orange', 5:'black'}
    for i_dim, dim in enumerate(dimension):
        fig = plt.figure(figsize=(7, 5))
        ax = fig.add_subplot(111)
        yranges = [np.inf, -np.inf]
        legend_handles = []
        xticks = [r'$\mathcal{O}_{12}$', r'$\mathcal{O}_{13}$', r'$\mathcal{O}_{23}$']
        ax.set_xlim(0, len(xticks)+1)
        ax.set_xticklabels([''] + xticks + [''])
        ax.set_xlabel(r'BSM operator angle')
        ylabel = r'${\rm log}_{10} \Lambda' + get_units(dim) + r'$'
        ax.set_ylabel(ylabel)
        for i_th in xrange(len(xticks)):
            for i_frs, frs in enumerate(fix_sfr_mfr):
                scale, evidences = angles_lim_array[i_dim][i_frs][i_th].T
                null = evidences[0]
                # TODO(shivesh): negative or not?
                reduced_ev = -(evidences - null)
                al = scale[reduced_ev > significance]
                if len(al) > 0:
                    label = '[{0}, {1}, {2}]'.format(frs[3], frs[4], frs[5])
                    lim = al[0]
                    print 'frs, dim, lim = ', frs, dim, lim
                    if lim < yranges[0]: yranges[0] = lim
                    if lim > yranges[1]: yranges[1] = lim+4
                    line = plt.Line2D(
                        (i_th+1-0.1, i_th+1+0.1), (lim, lim), lw=3, color=colour[i_frs], label=label
                    )
                    ax.add_line(line)
                    if i_th == 0: legend_handles.append(line)
                    x_offset = i_frs*0.05 - 0.05
                    ax.annotate(
                        s='', xy=(i_th+1+x_offset, lim), xytext=(i_th+1+x_offset, lim+3),
                        arrowprops={'arrowstyle': '<-', 'lw': 1.2, 'color':colour[i_frs]}
                    )
        try:
            yranges = (myround(yranges[0], up=True), myround(yranges[1], down=True))
            ax.set_ylim(yranges)
        except: pass

        ax.legend(handles=legend_handles, prop=dict(size=8), loc='upper right',
                  title='dimension {0}'.format(dim))

        for ymaj in ax.yaxis.get_majorticklocs():
            ax.axhline(y=ymaj, ls='-', color='gray', alpha=0.4, linewidth=1)
        for xmaj in ax.xaxis.get_majorticklocs():
            ax.axvline(x=xmaj, ls='-', color='gray', alpha=0.4, linewidth=1)

        for of in outformat:
            fig.savefig('./images/angles_limit_DIM{0}'.format(dim)+'.'+of, bbox_inches='tight', dpi=150)
