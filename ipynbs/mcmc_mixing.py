#! /usr/bin/env python
"""
Sample points from a gaussian likelihood
"""

from __future__ import absolute_import, division

import sys

import argparse
import multiprocessing
from fractions import gcd

import numpy as np
from scipy.stats import multivariate_normal

import emcee
import tqdm

DTYPE  = np.float128
CDTYPE = np.complex128


def solve_ratio(fr):
    denominator = reduce(gcd, fr)
    return [int(x/denominator) for x in fr]


def normalise_fr(fr):
    return np.array(fr) / float(np.sum(fr))

SOURCE_FR = normalise_fr([1, 1, 1])
SIGMA = 0.001
ANGLES = (np.sin(np.pi/4.)**2, 1.0, 0, 0)


def angles_to_u(bsm_angles):
    """Convert angular projection of the mixing matrix elements back into the
    mixing matrix elements.

    Parameters
    ----------
    bsm_angles : list, length = 4
        sin(12)^2, cos(13)^4, sin(23)^2 and deltacp

    Returns
    ----------
    unitary numpy ndarray of shape (3, 3)

    Examples
    ----------
    >>> from fr import angles_to_u
    >>> print angles_to_u((0.2, 0.3, 0.5, 1.5))
    array([[ 0.66195018+0.j        ,  0.33097509+0.j        ,  0.04757188-0.6708311j ],
           [-0.34631487-0.42427084j,  0.61741198-0.21213542j,  0.52331757+0.j        ],
           [ 0.28614067-0.42427084j, -0.64749908-0.21213542j,  0.52331757+0.j        ]])

    """
    s12_2, c13_4, s23_2, dcp = bsm_angles
    dcp = CDTYPE(dcp)

    c12_2 = 1. - s12_2
    c13_2 = np.sqrt(c13_4, dtype=DTYPE)
    s13_2 = 1. - c13_2
    c23_2 = 1. - s23_2

    t12 = np.arcsin(np.sqrt(s12_2, dtype=DTYPE))
    t13 = np.arccos(np.sqrt(c13_2, dtype=DTYPE))
    t23 = np.arcsin(np.sqrt(s23_2, dtype=DTYPE))

    c12 = np.cos(t12)
    s12 = np.sin(t12)
    c13 = np.cos(t13)
    s13 = np.sin(t13)
    c23 = np.cos(t23)
    s23 = np.sin(t23)

    p1 = np.array([[1   , 0   , 0]                   , [0    , c23 , s23] , [0                   , -s23 , c23]] , dtype=CDTYPE)
    p2 = np.array([[c13 , 0   , s13*np.exp(-1j*dcp)] , [0    , 1   , 0]   , [-s13*np.exp(1j*dcp) , 0    , c13]] , dtype=CDTYPE)
    p3 = np.array([[c12 , s12 , 0]                   , [-s12 , c12 , 0]   , [0                   , 0    , 1]]   , dtype=CDTYPE)

    u = np.dot(np.dot(p1, p2), p3)
    return u


def u_to_fr(source_fr, matrix):
    """Compute the observed flavour ratio assuming decoherence.

    Parameters
    ----------
    source_fr : list, length = 3
        Source flavour ratio components

    matrix : numpy ndarray, dimension 3
        Mixing matrix

    Returns
    ----------
    Measured flavour ratio

    Examples
    ----------
    >>> from fr import params_to_BSMu, u_to_fr
    >>> print u_to_fr((1, 2, 0), params_to_BSMu((0.2, 0.3, 0.5, 1.5, -20), 3, 1000))
        array([ 0.33740075,  0.33176584,  0.33083341])

    """
    composition = np.einsum(
        'ai, bi, a -> b', np.abs(matrix)**2, np.abs(matrix)**2, source_fr,
    )
    ratio = composition / np.sum(source_fr)
    return ratio


MIXING_U = angles_to_u(ANGLES)
MEASURED_FR = u_to_fr(SOURCE_FR, MIXING_U)


def angles_to_fr(angles):
    sphi4, c2psi = angles

    psi = (0.5)*np.arccos(c2psi)

    sphi2 = np.sqrt(sphi4)
    cphi2 = 1. - sphi2
    spsi2 = np.sin(psi)**2
    cspi2 = 1. - spsi2

    x = sphi2*cspi2
    y = sphi2*spsi2
    z = cphi2
    return x, y, z


def triangle_llh(theta):
    """-Log likelihood function for a given theta."""
    fr = angles_to_fr(theta)
    fr_bf = MEASURED_FR
    cov_fr = np.identity(3) * SIGMA
    return np.log(multivariate_normal.pdf(fr, mean=fr_bf, cov=cov_fr))


def lnprior(theta):
    """Priors on theta."""
    sphi4, c2psi = theta

    # Flavour ratio bounds
    if 0. <= sphi4 <= 1.0 and -1.0 <= c2psi <= 1.0:
        pass
    else: return -np.inf

    return 0.


def lnprob(theta):
    """Prob function for mcmc."""
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + triangle_llh(theta)


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--source-ratio', type=int, nargs=3, default=[1, 1, 1],
        help='Set the central value for the source flavour ratio at IceCube'
    )
    parser.add_argument(
        '--angles', type=float, nargs=3, default=[0.307, (1-0.02195)**2, 0.565],
        help='Set the 1 sigma for the measured flavour ratio'
    )
    parser.add_argument(
        '--sigma-ratio', type=float, default=0.01,
        help='Set the 1 sigma for the measured flavour ratio'
    )
    parser.add_argument(
        '--burnin', type=int, default=100,
        help='Amount to burnin'
    )
    parser.add_argument(
        '--nwalkers', type=int, default=100,
        help='Number of walkers'
    )
    parser.add_argument(
        '--nsteps', type=int, default=2000,
        help='Number of steps to run'
    )
    parser.add_argument(
        '--seed', type=int, default=99,
        help='Set the random seed value'
    )
    parser.add_argument(
        '--outfile', type=str, default='./untitled',
        help='Path to output chains'
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    np.random.seed(args.seed)

    global SOURCE_FR
    global SIGMA
    global ANGLES
    global MIXING_U
    global MEASURED_FR

    SOURCE_FR = np.array(args.source_ratio) / float(np.sum(args.source_ratio))
    SIGMA = args.sigma_ratio
    ANGLES = args.angles + [0]
    MIXING_U = angles_to_u(ANGLES)
    MEASURED_FR = u_to_fr(SOURCE_FR, MIXING_U)

    print 'SOURCE_FR = {0}'.format(SOURCE_FR)
    print 'SIGMA = {0}'.format(SIGMA)
    print 'ANGLES = {0}'.format(ANGLES)
    print 'MIXING_U = {0}'.format(MIXING_U)
    print 'MEASURED_FR = {0}'.format(MEASURED_FR)

    ndim = 2
    nwalkers = args.nwalkers
    ntemps = 1
    burnin = args.burnin
    betas = np.array([1e0, 1e-1, 1e-2, 1e-3, 1e-4])
    p0_base = [0.5, 0.5]
    p0_std = [.2, 0.2]

    print 'p0_base', p0_base
    print 'p0_std', p0_std
    p0 = np.random.normal(p0_base, p0_std, size=[ntemps, nwalkers, ndim])
    print map(lnprior, p0[0])

    # threads = multiprocessing.cpu_count()
    threads = 1
    sampler = emcee.PTSampler(
        ntemps, nwalkers, ndim, triangle_llh, lnprior, threads=threads
    )

    print "Running burn-in"
    for result in tqdm.tqdm(sampler.sample(p0, iterations=burnin), total=burnin):
        pos, prob, state = result
    sampler.reset()
    print "Finished burn-in"

    nsteps = args.nsteps

    print "Running"
    for _ in tqdm.tqdm(sampler.sample(pos, iterations=nsteps), total=nsteps):
        pass
    print "Finished"

    sr = solve_ratio(SOURCE_FR)
    outfile = args.outfile+'_{0}_{1}_{2}_{3:.1E}_{4:.2f}_{5:.2f}_{6:.2f}'.format(
            sr[0], sr[1], sr[2], SIGMA, ANGLES[0], ANGLES[1], ANGLES[2]
    )

    samples = sampler.chain[0, :, :, :].reshape((-1, ndim))
    print 'acceptance fraction', sampler.acceptance_fraction
    print 'sum of acceptance fraction', np.sum(sampler.acceptance_fraction)
    print 'np.unique(samples[:,0]).shape', np.unique(samples[:,0]).shape

    try:
        print 'autocorrelation', sampler.acor
    except:
        print 'WARNING : NEED TO RUN MORE SAMPLES FOR FILE ' + outfile
    print 'outfile = ', outfile

    fr_samples = np.array(map(angles_to_fr, samples))
    np.save(outfile+'.npy', fr_samples)

    print "DONE!"


main.__doc__ = __doc__


if __name__ == '__main__':
    main()

