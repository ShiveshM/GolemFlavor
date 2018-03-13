#! /usr/bin/env python
"""
From an MCMC chains file, make a triangle plot.
"""

from __future__ import absolute_import, division

import sys, os
import errno

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import rc, rcParams

import getdist
from getdist import plots
from getdist import mcsamples

import mcmc_scan

rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern'], 'size':18})
cols = ['#29A2C6','#FF6D31','#FFCB18','#73B66B','#EF597B', '#333333']

font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 18}


def plot(infile, angles, outfile, measured_ratio, sigma_ratio, fix_sfr,
         fix_mixing, fix_scale, source_ratio, scale, dimension, energy, scale_bounds):
    """Make the triangle plot"""
    if not angles:
        if fix_mixing:
            labels = [r'{\rm log}_{10}\Lambda', r'\phi_e', r'\phi_\mu', r'\phi_\tau']
        elif fix_sfr:
            if fix_scale:
                labels = [r'\mid \tilde{U}_{e1} \mid', r'\mid \tilde{U}_{e2} \mid', r'\mid \tilde{U}_{e3} \mid', \
                          r'\mid \tilde{U}_{\mu1} \mid', r'\mid \tilde{U}_{\mu2} \mid', r'\mid \tilde{U}_{\mu3} \mid', \
                          r'\mid \tilde{U}_{\tau1} \mid', r'\mid \tilde{U}_{\tau2} \mid', r'\mid \tilde{U}_{\tau3} \mid']
            else:
                labels = [r'\mid \tilde{U}_{e1} \mid', r'\mid \tilde{U}_{e2} \mid', r'\mid \tilde{U}_{e3} \mid', \
                          r'\mid \tilde{U}_{\mu1} \mid', r'\mid \tilde{U}_{\mu2} \mid', r'\mid \tilde{U}_{\mu3} \mid', \
                          r'\mid \tilde{U}_{\tau1} \mid', r'\mid \tilde{U}_{\tau2} \mid', r'\mid \tilde{U}_{\tau3} \mid', \
                          r'{\rm log}_{10}(\Lambda)']
        else:
            if fix_scale:
                labels = [r'\mid \tilde{U}_{e1} \mid', r'\mid \tilde{U}_{e2} \mid', r'\mid \tilde{U}_{e3} \mid', \
                          r'\mid \tilde{U}_{\mu1} \mid', r'\mid \tilde{U}_{\mu2} \mid', r'\mid \tilde{U}_{\mu3} \mid', \
                          r'\mid \tilde{U}_{\tau1} \mid', r'\mid \tilde{U}_{\tau2} \mid', r'\mid \tilde{U}_{\tau3} \mid', \
                          r'{\rm log}_{10}\Lambda', r'\phi_e', r'\phi_\mu', r'\phi_\tau']
            else:
                labels = [r'\mid \tilde{U}_{e1} \mid', r'\mid \tilde{U}_{e2} \mid', r'\mid \tilde{U}_{e3} \mid', \
                          r'\mid \tilde{U}_{\mu1} \mid', r'\mid \tilde{U}_{\mu2} \mid', r'\mid \tilde{U}_{\mu3} \mid', \
                          r'\mid \tilde{U}_{\tau1} \mid', r'\mid \tilde{U}_{\tau2} \mid', r'\mid \tilde{U}_{\tau3} \mid', \
                          r'\phi_e', r'\phi_\mu', r'\phi_\tau']
    else:
        if fix_sfr:
            if fix_mixing:
                assert 0
                labels=[r'\tilde{s}_{12}^2', r'{\rm log}_{10}\Lambda']
            elif fix_scale:
                labels=[r'\tilde{s}_{12}^2', r'\tilde{c}_{13}^4',
                        r'\tilde{s}_{23}^2', r'\tilde{\delta_{CP}}']
            else:
                labels=[r'\tilde{s}_{12}^2', r'\tilde{c}_{13}^4',
                        r'\tilde{s}_{23}^2', r'\tilde{\delta_{CP}}',
                        r'{\rm log}_{10}\Lambda']
        else:
            if fix_mixing:
                labels=[r'{\rm log}_{10}\Lambda', r'sin^4(\phi)', r'cos(2\psi)']
            elif fix_scale:
                labels=[r'\tilde{s}_{12}^2', r'\tilde{c}_{13}^4',
                        r'\tilde{s}_{23}^2', r'\tilde{\delta_{CP}}',
                        r'sin^4(\phi)', r'cos(2\psi)']
            else:
                labels=[r'\tilde{s}_{12}^2', r'\tilde{c}_{13}^4',
                        r'\tilde{s}_{23}^2', r'\tilde{\delta_{CP}}',
                        r'{\rm log}_{10}\Lambda', r'sin^4(\phi)', r'cos(2\psi)']
    labels = [r'convNorm', r'promptNorm', 'muonNorm', 'astroNorm', 'astroDeltaGamma'] + labels
    print 'labels', labels

    if not fix_scale:
        s2 = np.log10(scale_bounds)

    if not angles:
        if fix_mixing:
            ranges = [s2, (0, 1), (0, 1), (0, 1)]
        elif fix_sfr:
            if fix_scale:
                ranges = [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1)]
            else:
                ranges = [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), s2]
        else:
            if fix_scale:
                ranges = [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1)]
            else:
                ranges = [(0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), (0, 1), s2, (0, 1), (0, 1), (0, 1)]
    else:
        if fix_sfr:
            if fix_mixing:
                ranges = [(0, 1), s2]
            elif fix_scale:
                ranges = [(0, 1), (0, 1), (0, 1), (0, 2*np.pi)]
            else:
                ranges = [(0, 1), (0, 1), (0, 1), (0, 2*np.pi), s2]
        else:
            if fix_mixing:
                ranges = [s2, (0, 1), (-1, 1)]
            elif fix_scale:
                ranges = [(0, 1), (0, 1), (0, 1), (0, 2*np.pi), (0, 1), (-1, 1)]
            else:
                ranges = [(0, 1), (0, 1), (0, 1), (0, 2*np.pi), s2, (0, 1), (-1, 1)]
    ranges = [(0, 5), (0, 5), (0, 5), (0, 5), (0, 5)] + ranges
    print 'ranges', ranges

    def flat_angles_to_u(x):
        return abs(mcmc_scan.angles_to_u(x)).astype(np.float32).flatten().tolist()

    raw = np.load(infile)
    print 'raw.shape', raw.shape
    if not angles:
        nuisance, raw = raw[:,5:], raw[:,-5:]
        if fix_mixing:
            fr_elements = np.array(map(mcmc_scan.angles_to_fr, raw[:,-2:]))
            sc_elements = raw[:,:-2]
            Tchain = np.column_stack([sc_elements, fr_elements])
        elif fix_sfr:
            if fix_scale:
                Tchain = np.array(map(flat_angles_to_u, raw))
            else:
                sc_elements = raw[:,-1:]
                m_elements = np.array(map(flat_angles_to_u, raw[:,:-1]))
                Tchain = np.column_stack([m_elements, sc_elements])
        else:
            if fix_scale:
                fr_elements = np.array(map(mcmc_scan.angles_to_fr, raw[:,-2:]))
                m_elements = np.array(map(flat_angles_to_u, raw[:,:-2]))
                Tchain = np.column_stack([m_elements, fr_elements])
            else:
                fr_elements = np.array(map(mcmc_scan.angles_to_fr, raw[:,-2:]))
                sc_elements = raw[:,-3:-2]
                m_elements = np.array(map(flat_angles_to_u, raw[:,:-3]))
                Tchain = np.column_stack([m_elements, sc_elements, fr_elements])
        Tchain = np.column_stack([nuisance, Tchain])
    else:
        Tchain = raw
    print 'Tchain.shape', Tchain.shape

    if fix_sfr:
        if fix_scale:
            label = 'Source flavour ratio = [{0:.2f}, {1:.2f}, {2:.2f}]\nIC observed flavour ratio = [{3:.2f}, {4:.2f}, {5:.2f}]\nSigma = {6:.3f}\nDimension = {7}\nEnergy = {8} GeV\nScale = {9}'.format(
                source_ratio[0], source_ratio[1], source_ratio[2],
                measured_ratio[0], measured_ratio[1], measured_ratio[2], sigma_ratio,
                dimension, int(energy), scale
            )
        else:
            label = 'Source flavour ratio = [{0:.2f}, {1:.2f}, {2:.2f}]\nIC observed flavour ratio = [{3:.2f}, {4:.2f}, {5:.2f}]\nSigma = {6:.3f}\nDimension = {7}\nEnergy = {8} GeV'.format(
                source_ratio[0], source_ratio[1], source_ratio[2],
                measured_ratio[0], measured_ratio[1], measured_ratio[2], sigma_ratio,
                dimension, int(energy)
            )
    else:
        if fix_scale:
	    label = 'IC observed flavour ratio = [{0:.2f}, {1:.2f}, {2:.2f}]\nSigma = {3:.3f}\nDimension = {4}\nEnergy = {5} GeV\nScale = {6}'.format(
		measured_ratio[0], measured_ratio[1], measured_ratio[2], sigma_ratio,
		dimension, int(energy), scale
	    )
	else:
	    label = 'IC observed flavour ratio = [{0:.2f}, {1:.2f}, {2:.2f}]\nSigma = {3:.3f}\nDimension = {4}\nEnergy = {5} GeV'.format(
		measured_ratio[0], measured_ratio[1], measured_ratio[2], sigma_ratio,
		dimension, int(energy)
	    )

    Tsample = mcsamples.MCSamples(
        samples=Tchain, labels=labels, ranges=ranges
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
    if fix_mixing and fix_sfr:
        mpl.pyplot.figtext(0.4, 0.7, label, fontsize=4)
    else:
        mpl.pyplot.figtext(0.5, 0.7, label, fontsize=15)
    print 'outfile = {0}'.format(outfile)
    try:
        os.makedirs(outfile[:-len(os.path.basename(outfile))])
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(outfile[:-len(os.path.basename(outfile))]):
            pass
        else:
            raise
    g.export(outfile)


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--infile', type=str, required=True,
        help='Path to MCMC chains'
    )
    parser.add_argument(
        '--angles', default=False, action='store_true',
        help='Plot in terms of mixing angles'
    )
    parser.add_argument(
        '--outfile', type=str, default='./untitled.pdf',
        help='Path to output plot'
    )
    parser.add_argument(
        '--bestfit-ratio', type=int, nargs=3, required=False,
        help='Set the bestfit flavour ratio'
    )
    parser.add_argument(
        '--sigma-ratio', type=float, required=False,
        help='Set the 1 sigma for the flavour ratio'
    )
    parser.add_argument(
        '--fix-sfr', action='store_true',
        help='Fix the source flavour ratio'
    )
    parser.add_argument(
        '--fix-mixing', action='store_true',
        help='Fix the new physics mixing values to a single term, s_12^2'
    )
    parser.add_argument(
        '--source-ratio', type=int, nargs=3, default=[2, 1, 0],
        help='Set the source flavour ratio for the case when you want to fix it'
    )
    parser.add_argument(
        '--scale', type=float, required=False,
        help='Fix the scale to this value'
    )
    parser.add_argument(
        '--dimension', type=int, default=3, help='Dimension'
    )
    parser.add_argument(
        '--energy', type=float, default=1000, help='Energy'
    )
    parser.add_argument(
        '--scale-bounds', type=float, nargs=2,
        help='Upper and lower limits to plot the new physics scale'
    )
    args = parser.parse_args()
    return args


def main():
    args = vars(parse_args())
    plot(**args)

    print "DONE!"


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
