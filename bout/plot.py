import os

import numpy as np
import numpy.ma as ma

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.offsetbox import AnchoredText
from matplotlib import rc

rc('text', usetex=False)
rc('font', **{'family':'serif', 'serif':['Computer Modern'], 'size':18})

fix_sfr_mfr = [
    (1, 1, 1, 1, 2, 0),
    # (1, 1, 1, 1, 0, 0),
    (1, 1, 1, 0, 1, 0),
]

# FR
# dimension         = [3, 6]
dimension         = [3, 6]
sigma_ratio       = ['0.01']
energy_dependance = 'spectral'
spectral_index    = -2
binning           = [1e4, 1e7, 5]
fix_mixing        = 'False'
fix_mixing_almost = 'False'
scale_region      = "1E10"

# Likelihood
likelihood = 'golemfit'
confidence = 2.71 # 90% for 1DOF
outformat = ['png']


def gen_identifier(measured_ratio, source_ratio, dimension, sigma_ratio=0.01):
    mr = np.array(measured_ratio) / float(np.sum(measured_ratio))
    sr = np.array(source_ratio) / float(np.sum(source_ratio))
    si = sigma_ratio
    out = '_{0:03d}_{1:03d}_{2:03d}_{3:04d}_sfr_{4:03d}_{5:03d}_{6:03d}_DIM{7}_single_scale'.format(
        int(mr[0]*100), int(mr[1]*100), int(mr[2]*100), int(si*1000),
        int(sr[0]*100), int(sr[1]*100), int(sr[2]*100), dimension
    )
    return out


def get_units(dimension):
    if dimension == 3: return r' / GeV'
    if dimension == 4: return r''
    if dimension == 5: return r' / GeV^{-1}'
    if dimension == 6: return r' / GeV^{-2}'
    if dimension == 7: return r' / GeV^{-3}'
    if dimension == 8: return r' / GeV^{-4}'


def myround(x, base=5, up=False, down=False):
    if up == down and up is True: assert 0
    if up: return int(base * np.round(float(x)/base-0.5))
    elif down: return int(base * np.round(float(x)/base+0.5))
    else: int(base * np.round(float(x)/base))


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
    for i_frs, frs in enumerate(fix_sfr_mfr):
        print '== DIM{0}'.format(dim)
        print '== FRS = {0}'.format(frs)
        outchain_head = '/data/user/smandalia/flavour_ratio/data/{0}/DIM{1}/SI_{2}/fix_ifr/0_01/'.format(likelihood, dim, spectral_index)
        infile = outchain_head + '/angles_limit/fr_anfr_evidence'+ gen_identifier(frs[:3], frs[-3:], dim) + '.npy'
        try:
            array = np.load(infile)
        except IOError:
            print 'failed to open {0}'.format(infile)
            continue
        print 'array', array
        print 'array', array.shape
        for i_th in xrange(len(xticks)):
            scale, llhs = array[i_th].T
            min_llh = np.min(llhs)
            delta_llh = 2*(llhs - min_llh)
            print 'scale', scale
            print 'delta_llh', delta_llh
            al = scale[delta_llh < confidence]
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
            else:
                print 'No points for DIM {0} FRS {1} NULL {2}!'.format(dim, frs, min_llh)
    try:
        yranges = (myround(yranges[0], up=True), myround(yranges[1], down=True))
        # ax.set_ylim(yranges)
        ax.set_ylim([-30, -20])
    except: pass

    ax.legend(handles=legend_handles, prop=dict(size=8), loc='upper right',
              title='dimension {0}'.format(dim))
    for ymaj in ax.yaxis.get_majorticklocs():
        ax.axhline(y=ymaj, ls=':', color='gray', alpha=0.4, linewidth=1)
    for xmaj in ax.xaxis.get_majorticklocs():
        ax.axvline(x=xmaj, ls=':', color='gray', alpha=0.4, linewidth=1)

    for of in outformat:
        fig.savefig('../images/freq/lim_DIM{0}.'.format(dim)+of, bbox_inches='tight', dpi=150)
