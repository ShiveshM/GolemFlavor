# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : March 17, 2018

"""
Useful functions to use an MCMC for the BSM flavour ratio analysis
"""

from __future__ import absolute_import, division

import argparse
from functools import partial

import emcee
import tqdm

import numpy as np
from scipy.stats import multivariate_normal

from utils import fr as fr_utils
from utils.enums import Likelihood
from utils.misc import enum_keys, enum_parse, make_dir, parse_bool


def lnprior(theta, paramset):
    """Priors on theta."""
    ranges = paramset.ranges
    for value, range in zip(theta, ranges):
        if range[0] <= value <= range[1]:
            pass
        else: return -np.inf
    return 0.


def mcmc(p0, triangle_llh, lnprior, ndim, nwalkers, burnin, nsteps, ntemps=1, threads=1):
    """Run the MCMC."""
    sampler = emcee.PTSampler(
        ntemps, nwalkers, ndim, triangle_llh, lnprior, threads=threads
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

    samples = sampler.chain[0, :, :, :].reshape((-1, ndim))
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
        '--likelihood', default='gaussian', type=partial(enum_parse, c=Likelihood),
        choices=enum_keys(Likelihood), help='likelihood contour'
    )
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


def gaussian_llh(fr, fr_bf, sigma):
    """Multivariate gaussian likelihood."""
    cov_fr = np.identity(3) * sigma
    return np.log(multivariate_normal.pdf(fr, mean=fr_bf, cov=cov_fr))


def gaussian_seed(paramset, ntemps, nwalkers):
    """Get gaussian seed values for the MCMC."""
    ndim = len(paramset)
    p0 = np.random.normal(
        paramset.values, paramset.stds, size=[ntemps, nwalkers, ndim]
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


def triangle_llh(theta, args):
    """-Log likelihood function for a given theta."""
    if args.fix_source_ratio:
        fr1, fr2, fr3 = args.source_ratio
        bsm_angles = theta[-5:]
    else:
        fr1, fr2, fr3 = fr_utils.angles_to_fr(theta[-2:])
        bsm_angles = theta[-7:-2]

    u = fr_utils.params_to_BSMu(
        theta      = bsm_angles,
        dim        = args.dimension,
        energy     = args.energy,
        no_bsm     = args.no_bsm,
        fix_mixing = args.fix_mixing,
        fix_scale  = args.fix_scale,
        scale      = args.scale
    )

    fr = fr_utils.u_to_fr((fr1, fr2, fr3), u)
    fr_bf = args.measured_ratio
    if args.likelihood is Likelihood.FLAT:
        return 1.
    elif args.likelihood is Likelihood.GAUSSIAN:
        return gaussian_llh(fr, fr_bf, args.sigma_ratio)
    elif args.likelihood is Likelihood.GOLEMFIT:
        raise NotImplementedError('TODO')
        import GolemFitPy as gf
        from collections import namedtuple
        datapaths = gf.DataPaths()
        IceModels = namedtuple('IceModels', 'std_half2')
        fitters = IceModels(*[
            gf.GolemFit(datapaths,
                        gf.gen_steering_params(SteeringCateg.__members__[ice]),
                        xs_params(XSCateg.nom)) for ice in IceModels._fields])
        for fitter in fitters:
            fitter.SetFitParametersFlag(fit_flags(FitFlagCateg.xs))
            fitter.SetFitPriors(fit_priors(FitPriorsCateg.xs))

