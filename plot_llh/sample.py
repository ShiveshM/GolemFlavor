#! /usr/bin/env python
"""
Sample points for a specific scenario
"""

from __future__ import absolute_import, division

import sys
sys.path.extend(['.', '../'])

import argparse
from functools import partial

import numpy as np
from scipy.stats import multivariate_normal

from utils import fr as fr_utils
from utils import likelihood as llh_utils
from utils import mcmc as mcmc_utils
from utils import misc as misc_utils
from utils.param import Param, ParamSet, get_paramsets
from utils.enums import MixingScenario, MCMCSeedType


def triangle_llh(theta, args):
    """-Log likelihood function for a given theta."""
    fr = fr_utils.angles_to_fr(theta)
    fr_bf = args.measured_ratio
    cov_fr = np.identity(3) * args.sigma_ratio
    return np.log(multivariate_normal.pdf(fr, mean=fr_bf, cov=cov_fr))


def ln_prior(theta):
    """Priors on theta."""
    sphi4, c2psi = theta
    # Flavour ratio bounds
    if 0. <= sphi4 <= 1.0 and -1.0 <= c2psi <= 1.0:
        pass
    else: return -np.inf
    return 0.


def ln_prob(theta, args):
    """Prob function for mcmc."""
    lp = ln_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + triangle_llh(theta, args)


def process_args(args):
    """Process the input args."""
    if args.fix_mixing is MixingScenario.NONE:
        raise AssertionError('Must set a mixing scenario using --fix-mixing')
    if not args.fix_source_ratio:
        raise AssertionError('Must set source ratio using --fix-source-ratio')

    args.source_ratio = fr_utils.normalise_fr(args.source_ratio)
    if args.fix_mixing is MixingScenario.T12:
        s12_2, c13_4, s23_2, dcp = 0.5, 1.0, 0.0, 0.
    elif args.fix_mixing is MixingScenario.T13:
        s12_2, c13_4, s23_2, dcp = 0.0, 0.25, 0.0, 0.
    elif args.fix_mixing is MixingScenario.T23:
        s12_2, c13_4, s23_2, dcp = 0.0, 1.0, 0.5, 0.

    mm = np.array(
        fr_utils.angles_to_u((s12_2, c13_4, s23_2, dcp)), dtype=np.complex256
    )
    args.measured_ratio = fr_utils.u_to_fr(args.source_ratio, mm)

    if not args.fix_scale:
        args.scale, args.scale_region = fr_utils.estimate_scale(args)


def parse_args(args=None):
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="BSM flavour ratio analysis",
        formatter_class=misc_utils.SortingHelpFormatter,
    )
    parser.add_argument(
        '--seed', type=misc_utils.seed_parse, default='25',
        help='Set the random seed value'
    )
    parser.add_argument(
        '--threads', type=misc_utils.thread_type, default='1',
        help='Set the number of threads to use (int or "max")'
    )
    parser.add_argument(
        '--outfile', type=str, default='./untitled',
        help='Path to output chains'
    )
    parser.add_argument(
        '--plot-statistic', type=misc_utils.parse_bool, default='False',
        help='Plot MultiNest evidence or LLH value'
    )
    fr_utils.fr_argparse(parser)
    llh_utils.likelihood_argparse(parser)
    mcmc_utils.mcmc_argparse(parser)
    if args is None: return parser.parse_args()
    else: return parser.parse_args(args.split())


def main():
    args = parse_args()
    process_args(args)
    misc_utils.print_args(args)

    if args.seed is not None:
        np.random.seed(args.seed)

    paramset = ParamSet([
        Param(name='sphi4', value=0.5, ranges=[0.0,  1.0], std=0.2),
        Param(name='c2psi', value=0.0, ranges=[-1.0, 1.0], std=0.2)
    ])

    outfile = misc_utils.gen_outfile_name(args)
    print '== {0:<25} = {1}'.format('outfile', outfile)

    if args.run_mcmc:
        ndim = len(paramset)
        print paramset

        if args.mcmc_seed_type == MCMCSeedType.UNIFORM:
            p0 = mcmc_utils.flat_seed(paramset, nwalkers=args.nwalkers)

        samples = mcmc_utils.mcmc(
            p0       = p0,
            ln_prob  = partial(ln_prob, args=args),
            ndim     = ndim,
            nwalkers = args.nwalkers,
            burnin   = args.burnin,
            nsteps   = args.nsteps,
            threads  = args.threads
        )
        fr_chains = np.array(map(fr_utils.angles_to_fr, samples), dtype=float)
        mcmc_utils.save_chains(fr_chains, outfile)
    print "DONE!"


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
