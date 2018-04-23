# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : March 17, 2018

"""
Useful functions to use an MCMC for the BSM flavour ratio analysis
"""

from __future__ import absolute_import, division

from functools import partial

import emcee
import tqdm

import numpy as np

from utils.enums import MCMCSeedType
from utils.misc import enum_parse, make_dir, parse_bool


def mcmc(p0, ln_prob, ndim, nwalkers, burnin, nsteps, threads=1):
    """Run the MCMC."""
    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, ln_prob, threads=threads
    )

    print "Running burn-in"
    for result in tqdm.tqdm(sampler.sample(p0, iterations=burnin), total=burnin):
        pos, prob, state = result
    sampler.reset()
    print "Finished burn-in"

    print "Running"
    for _ in tqdm.tqdm(sampler.sample(pos, iterations=nsteps), total=nsteps):
        pass
    print "Finished"

    samples = sampler.chain.reshape((-1, ndim))
    print 'acceptance fraction', sampler.acceptance_fraction
    print 'sum of acceptance fraction', np.sum(sampler.acceptance_fraction)
    print 'np.unique(samples[:,0]).shape', np.unique(samples[:,0]).shape
    try:
        print 'autocorrelation', sampler.acor
    except:
        print 'WARNING : NEED TO RUN MORE SAMPLES'

    return samples


def mcmc_argparse(parser):
    parser.add_argument(
        '--run-mcmc', type=parse_bool, default='True',
        help='Run the MCMC'
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
        '--mcmc-seed-type', default='uniform',
        type=partial(enum_parse, c=MCMCSeedType), choices=MCMCSeedType,
        help='Type of distrbution to make the initial MCMC seed'
    )
    parser.add_argument(
        '--plot-angles', type=misc_utils.parse_bool, default='True',
        help='Plot MCMC triangle in the angles space'
    )
    parser.add_argument(
        '--plot-elements', type=misc_utils.parse_bool, default='False',
        help='Plot MCMC triangle in the mixing elements space'
    )


def flat_seed(paramset, nwalkers):
    """Get gaussian seed values for the MCMC."""
    ndim = len(paramset)
    low = np.array(paramset.seeds).T[0]
    high = np.array(paramset.seeds).T[1]
    p0 = np.random.uniform(
        low=low, high=high, size=[nwalkers, ndim]
    )
    return p0


def gaussian_seed(paramset, nwalkers):
    """Get gaussian seed values for the MCMC."""
    ndim = len(paramset)
    p0 = np.random.normal(
        paramset.values, paramset.stds, size=[nwalkers, ndim]
    )
    return p0


def save_chains(chains, outfile):
    """Save the chains.

    Parameters
    ----------
    chains : numpy ndarray
        MCMC chains to save

    outfile : str
        Output file location of chains

    """
    make_dir(outfile)
    print 'Saving chains to location {0}'.format(outfile+'.npy')
    np.save(outfile+'.npy', chains)

