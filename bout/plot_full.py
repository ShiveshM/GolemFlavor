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
        outchain_head = '/data/user/smandalia/flavour_ratio/data/{0}/DIM{1}/SI_{2}/fix_ifr/0_01/'.format(likelihood, dim, spectral_index)
        infile = outchain_head + '/bayes_factor/fr_fr_evidence' + gen_identifier(frs[:3], frs[-3:], dim) + '.npy'
        try:
            array = np.load(infile)
        except IOError:
            print 'failed to open {0}'.format(infile)
            continue
        print 'array', array
        print 'array', array.shape
        scale, llhs = array.T
        print 'scale min', scale[np.argmin(llhs)]
        null = llhs[np.argmin(llhs)]
        # null = llhs[0]
        # TODO(shivesh): negative or not?
        reduced_ev = 2*(llhs - null)
        print 'reduced_ev', reduced_ev
        al = scale[reduced_ev < confidence]
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
            # print 'scales, reduced_ev', np.dstack([scale.data, reduced_ev.data])
yranges = (myround(yranges[0], up=True), myround(yranges[1], down=True))
ax.set_ylim(yranges)

ax.legend(handles=legend_handles, prop=dict(size=8), loc='upper right')
for ymaj in ax.yaxis.get_majorticklocs():
    ax.axhline(y=ymaj, ls=':', color='gray', alpha=0.4, linewidth=1)
for xmaj in ax.xaxis.get_majorticklocs():
    ax.axvline(x=xmaj, ls=':', color='gray', alpha=0.4, linewidth=1)

for of in outformat:
    fig.savefig('../images/freq/full_corr.'+of, bbox_inches='tight', dpi=150)
