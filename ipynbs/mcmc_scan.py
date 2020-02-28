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

MEASURED_FR = [1, 1, 1]
SIGMA = 0.001


def solve_ratio(fr):
    denominator = reduce(gcd, fr)
    return [int(x/denominator) for x in fr]


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
        '--measured-ratio', type=int, nargs=3, default=[1, 1, 1],
        help='Set the central value for the measured flavour ratio at IceCube'
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

    global MEASURED_FR
    global SIGMA

    MEASURED_FR = np.array(args.measured_ratio) / float(np.sum(args.measured_ratio))
    SIGMA = args.sigma_ratio

    print 'MEASURED_FR = {0}'.format(MEASURED_FR)
    print 'SIGMA = {0}'.format(SIGMA)

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

    mr = solve_ratio(MEASURED_FR)
    outfile = args.outfile+'_{0}_{1}_{2}_{3:.1E}'.format(
            mr[0], mr[1], mr[2], SIGMA
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

