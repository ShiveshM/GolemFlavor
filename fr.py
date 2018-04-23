#! /usr/bin/env python
# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : March 17, 2018

"""
HESE BSM flavour ratio MCMC analysis script
"""

from __future__ import absolute_import, division

import argparse
from functools import partial

import numpy as np

from utils import fr as fr_utils
from utils import gf as gf_utils
from utils import likelihood as llh_utils
from utils import mcmc as mcmc_utils
from utils import misc as misc_utils
from utils import plot as plot_utils
from utils.enums import EnergyDependance, Likelihood, MCMCSeedType, ParamTag
from utils.fr import estimate_scale, normalise_fr
from utils.param import Param, ParamSet, get_paramsets


def define_nuisance():
    """Define the nuisance parameters to scan over with default values,
    ranges and sigma.
    """
    tag = ParamTag.NUISANCE
    nuisance = ParamSet(
        Param(name='convNorm',        value=1.,  seed=[0.5, 2. ], ranges=[0. , 50.], std=0.3,  tag=tag),
        Param(name='promptNorm',      value=0.,  seed=[0. , 6. ], ranges=[0. , 50.], std=0.05, tag=tag),
        Param(name='muonNorm',        value=1.,  seed=[0.1, 2. ], ranges=[0. , 50.], std=0.1,  tag=tag),
        Param(name='astroNorm',       value=6.9, seed=[0.1, 10.], ranges=[0. , 50.], std=0.1,  tag=tag),
        Param(name='astroDeltaGamma', value=2.5, seed=[1. , 3. ], ranges=[-5., 5. ], std=0.1,  tag=tag)
    )
    return nuisance


def nuisance_argparse(parser):
    nuisance_paramset = define_nuisance()
    for parm in nuisance_paramset:
        parser.add_argument(
            '--'+parm.name, type=float, default=parm.value,
            help=parm.name+' to inject'
        )


def process_args(args):
    """Process the input args."""
    if args.fix_mixing and args.fix_scale:
        raise NotImplementedError('Fixed mixing and scale not implemented')
    if args.fix_mixing and args.fix_mixing_almost:
        raise NotImplementedError(
            '--fix-mixing and --fix-mixing-almost cannot be used together'
        )

    args.measured_ratio = normalise_fr(args.measured_ratio)
    if args.fix_source_ratio:
        args.source_ratio = normalise_fr(args.source_ratio)

    if args.energy_dependance is EnergyDependance.SPECTRAL:
        args.binning = np.logspace(
            np.log10(args.binning[0]), np.log10(args.binning[1]), args.binning[2]+1
        )

    if not args.fix_scale:
        args.scale = estimate_scale(args)
        args.scale_region = (args.scale/args.scale_region, args.scale*args.scale_region)


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
        help='Path to output results'
    )
    fr_utils.fr_argparse(parser)
    gf_utils.gf_argparse(parser)
    llh_utils.likelihood_argparse(parser)
    mcmc_utils.mcmc_argparse(parser)
    nuisance_argparse(parser)
    if args is None: return parser.parse_args()
    else: return parser.parse_args(args.split())


def main():
    args = parse_args()
    process_args(args)
    misc_utils.print_args(args)

    if args.seed is not None:
        np.random.seed(args.seed)

    asimov_paramset, llh_paramset = get_paramsets(args)
    outfile = misc_utils.gen_outfile_name(args)
    print '== {0:<25} = {1}'.format('outfile', outfile)

    if args.run_mcmc:
        if args.likelihood is Likelihood.GOLEMFIT:
            fitter = gf_utils.setup_fitter(args, asimov_paramset)
        else: fitter = None

        ln_prob = partial(
            llh_utils.ln_prob, args=args, fitter=fitter,
            asimov_paramset=asimov_paramset, llh_paramset=llh_paramset
        )

        ndim = len(llh_paramset)
        if args.mcmc_seed_type == MCMCSeedType.UNIFORM:
            p0 = mcmc_utils.flat_seed(
                llh_paramset, nwalkers=args.nwalkers
            )
        elif args.mcmc_seed_type == MCMCSeedType.GAUSSIAN:
            p0 = mcmc_utils.gaussian_seed(
                llh_paramset, nwalkers=args.nwalkers
            )

        samples = mcmc_utils.mcmc(
            p0       = p0,
            ln_prob  = ln_prob,
            ndim     = ndim,
            nwalkers = args.nwalkers,
            burnin   = args.burnin,
            nsteps   = args.nsteps,
            threads  = 1
            # TODO(shivesh): broken because you cannot pickle a GolemFitPy object
            # threads      = misc_utils.thread_factors(args.threads)[0]
        )
        mcmc_utils.save_chains(samples, outfile)

    plot_utils.chainer_plot(
        infile       = outfile+'.npy',
        outfile      = outfile[:5]+outfile[5:].replace('data', 'plots'),
        outformat    = ['pdf'],
        args         = args,
        llh_paramset = llh_paramset
    )
    print "DONE!"


main.__doc__ = __doc__


if __name__ == '__main__':
    main()

