import os

import numpy as np
import numpy.ma as ma

import scipy.interpolate as interpolate

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
confidence = 4.61 # 90% for 2DOF
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

labels = [r'$sin^2(2\mathcal{O}_{12})$', r'$sin^2(2\mathcal{O}_{13})$', r'$sin^2(2\mathcal{O}_{23})$']
for i_dim, dim in enumerate(dimension):
    for i_frs, frs in enumerate(fix_sfr_mfr):
        print '== DIM{0}'.format(dim)
        print '== FRS = {0}'.format(frs)
        outchain_head = '/data/user/smandalia/flavour_ratio/data/{0}/DIM{1}/SI_{2}/fix_ifr/0_01/'.format(likelihood, dim, spectral_index)
        infile = outchain_head + '/angles_corr/fr_co_evidence'+ gen_identifier(frs[:3], frs[-3:], dim) + '.npy'
        # infile = '../mnrun/fr_co_evidence_033_033_033_0010_sfr_033_066_000_DIM6_single_scale.npy'
        try:
            array = ma.masked_invalid(np.load(infile))
        except IOError:
            print 'failed to open {0}'.format(infile)
            continue
        print 'array', array
        print 'array', array.shape
        for i_scen in xrange(len(labels)):
            d = array[i_scen].reshape(len(array[i_scen])**2, 3)
            fig = plt.figure(figsize=(7, 5))
            ax = fig.add_subplot(111)
            xranges = [np.inf, -np.inf]
            legend_handles = []
            ax.set_ylim(0, 1)
            ax.set_ylabel(labels[i_scen])
            xlabel = r'${\rm log}_{10} \Lambda' + get_units(dim) + r'$'
            ax.set_xlabel(xlabel)

            x = d[:,0]
            y = d[:,1]
            z = d[:,2]

            print 'x', x
            print 'y', y
            print 'z', z
            null_idx = np.argmin(z)
            null = z[null_idx]
            print 'null = {0}, for scale = {1}'.format(null, x[null_idx])
            z = 2*(z - null)
            print 'scale', x
            print 'delta_llh', z

            # x_ = np.linspace(np.min(x), np.max(x), 30)
            # y_ = np.linspace(np.min(y), np.max(y), 30)
            # z_ = interpolate.gridddata((x, y), z, (x_[None,:], y_[:,None]), method='nearest')

            data = np.array([x, y, z, np.ones(x.shape)]).T
            print 'data', data
            data_clean = []
            for d in data:
                if not np.any(np.isnan(d)): data_clean.append(d)
            data = np.vstack(data_clean)
            sort_column = 3
            data_sorted = data[data[:,sort_column].argsort()]
            uni, c = np.unique(data[:,sort_column], return_counts=True)
            print uni, c
            print len(uni)
            print np.unique(c)
            
            n = len(uni)
            assert len(np.unique(c)) == 1
            c = c[0]
            col_array = []
            for col in data_sorted.T:
                col_array.append(col.reshape(n, c))
            col_array = np.stack(col_array)
            sep_arrays = []
            for x_i in xrange(n):
                sep_arrays.append(col_array[:,x_i])
            
            print len(sep_arrays)
            sep_arrays = sep_arrays[0][:3]
            print sep_arrays

            llh_90_percent = (sep_arrays[2] < confidence)
            data_90_percent = sep_arrays.T[llh_90_percent].T
            print 'data_90_percent', data_90_percent

            ax.tick_params(axis='x', labelsize=11)
            ax.tick_params(axis='y', labelsize=11)

            mini, maxi = np.min(x), np.max(x)
            ax.set_xlim((mini, maxi))
            ax.set_ylim(0, 1)
            ax.grid(b=False)
            
            x_v = data_90_percent[0].round(decimals=4)
            y_v = data_90_percent[1].round(decimals=4)
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
            
            for of in outformat:
                plt.savefig('../images/freq/lim_corr_DIM{0}_AN{1}_FRS{2}'.format(dim, i_scen, i_frs)+of, bbox_inches='tight', dpi=150)

