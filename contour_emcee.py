#! /usr/bin/env python
# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : November 26, 2018

"""
HESE flavour ratio contour
"""

from __future__ import absolute_import, division

import os
import argparse
from functools import partial

import numpy as np

from utils import fr as fr_utils
from utils import gf as gf_utils
from utils import llh as llh_utils
from utils import misc as misc_utils
from utils import mcmc as mcmc_utils
from utils import plot as plot_utils
from utils.enums import str_enum
from utils.enums import DataType, Likelihood, MCMCSeedType, ParamTag, PriorsCateg
from utils.param import Param, ParamSet, get_paramsets

from pymultinest import Analyzer, run


def define_nuisance():
    """Define the nuisance parameters."""
    nuisance = []
    tag = ParamTag.NUISANCE
    lg_prior = PriorsCateg.LIMITEDGAUSS
    nuisance.extend([
        # Param(name='convNorm',                  value=1.,   seed=[0.5,  2. ],  ranges=[0.1,   10.],   std=0.4, prior=lg_prior, tag=tag),
        # Param(name='promptNorm',                value=0.,   seed=[0.,   6. ],  ranges=[0.,    20.],   std=2.4, prior=lg_prior, tag=tag),
        Param(name='convNorm',                  value=1.,   seed=[0.5,  2. ],  ranges=[0.1,   10.],   std=0.4, tag=tag),
        Param(name='promptNorm',                value=0.,   seed=[0.,   6. ],  ranges=[0.,    20.],   std=2.4, tag=tag),
        Param(name='muonNorm',                  value=1.,   seed=[0.1,  2. ],  ranges=[0.,    10.],   std=0.1, tag=tag),
        Param(name='astroNorm',                 value=6.9,  seed=[0.,   5. ],  ranges=[0.,    20.],   std=1.5, tag=tag),
        Param(name='astroDeltaGamma',           value=2.5,  seed=[2.4,  3. ],  ranges=[-5.,   5. ],   std=0.1, tag=tag),
        Param(name='CRDeltaGamma',              value=0.,   seed=[-0.1, 0.1 ], ranges=[-1.,   1. ],   std=0.1, tag=tag),
        Param(name='NeutrinoAntineutrinoRatio', value=1.,   seed=[0.8,  1.2 ], ranges=[0.,    2. ],   std=0.1, tag=tag),
        Param(name='anisotropyScale',           value=1.,   seed=[0.8,  1.2 ], ranges=[0.,    2. ],   std=0.1, tag=tag),
        Param(name='domEfficiency',             value=0.99, seed=[0.8,  1.2 ], ranges=[0.8,   1.2 ],  std=0.1, tag=tag),
        Param(name='holeiceForward',            value=0.,   seed=[-0.8, 0.8 ], ranges=[-4.42, 1.58 ], std=0.1, tag=tag),
        Param(name='piKRatio',                  value=1.0,  seed=[0.8,  1.2 ], ranges=[0.,    2. ],   std=0.1, tag=tag)
    ])
    return ParamSet(nuisance)


def nuisance_argparse(parser):
    nuisance = define_nuisance()
    for parm in nuisance:
        parser.add_argument(
            '--'+parm.name, type=float, default=parm.value,
            help=parm.name+' to inject'
        )

def process_args(args):
    """Process the input args."""
    if args.likelihood is not Likelihood.GOLEMFIT \
       and args.likelihood is not Likelihood.GF_FREQ:
        raise AssertionError(
            'Likelihood method {0} not supported for this '
            'script!\nChoose either GOLEMFIT or GF_FREQ'.format(
                str_enum(args.likelihood)
            )
        )


def parse_args(args=None):
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="BSM flavour ratio analysis",
        formatter_class=misc_utils.SortingHelpFormatter,
    )
    parser.add_argument(
        '--injected-ratio', type=float, nargs=3, default=[1, 1, 1],
        help='Set the central value for the injected flavour ratio at IceCube'
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
    try:
        gf_utils.gf_argparse(parser)
    except: pass
    llh_utils.likelihood_argparse(parser)
    mcmc_utils.mcmc_argparse(parser)
    nuisance_argparse(parser)
    misc_utils.remove_option(parser, 'sigma_ratio')
    if args is None: return parser.parse_args()
    else: return parser.parse_args(args.split())


def gen_identifier(args):
    f = '_{0}_{1}'.format(*map(str_enum, (args.likelihood, args.data)))
    if args.data is not DataType.REAL:
        ir1, ir2, ir3 = misc_utils.solve_ratio(args.injected_ratio)
        f += '_INJ_{0:03d}_{1:03d}_{2:03d}'.format(ir1, ir2, ir3)
    return f


def gen_figtext(args, asimov_paramset):
    f = ''
    if args.data is DataType.REAL:
        f += 'IceCube Preliminary'
    else:
        ir1, ir2, ir3 = misc_utils.solve_ratio(args.injected_ratio)
        f += 'Injected ratio = [{0}, {1}, {2}]'.format(ir1, ir2, ir3)
        for param in asimov_paramset:
            f += '\nInjected {0:20s} = {1:.3f}'.format(
                param.name, param.nominal_value
            )
    return f


def triangle_llh(theta, args, hypo_paramset, fitter):
    """Log likelihood function for a given theta."""
    if len(theta) != len(hypo_paramset):
        raise AssertionError(
            'Dimensions of scan is not the same as the input '
            'params\ntheta={0}\nparamset]{1}'.format(theta, hypo_paramset)
        )
    for idx, param in enumerate(hypo_paramset):
        param.value = theta[idx]

    if args.likelihood is Likelihood.GOLEMFIT:
        llh = gf_utils.get_llh(fitter, hypo_paramset)
    elif args.likelihood is Likelihood.GF_FREQ:
        llh = gf_utils.get_llh_freq(fitter, hypo_paramset)

    return llh


def ln_prob(theta, args, hypo_paramset, fitter):
    lp = llh_utils.lnprior(theta, paramset=hypo_paramset)
    if not np.isfinite(lp):
        return -np.inf
    return lp + triangle_llh(
        theta,
        args           = args,
        hypo_paramset  = hypo_paramset,
        fitter         = fitter
    )


def main():
    args = parse_args()
    process_args(args)
    misc_utils.print_args(args)

    if args.seed is not None:
        np.random.seed(args.seed)

    asimov_paramset, hypo_paramset = get_paramsets(args, define_nuisance())
    hypo_paramset.extend(asimov_paramset.from_tag(ParamTag.BESTFIT))
    outfile = args.outfile + gen_identifier(args)
    print '== {0:<25} = {1}'.format('outfile', outfile)

    n_params = len(hypo_paramset)
    outfile = outfile + '_emcee_'

    print 'asimov_paramset', asimov_paramset
    print 'hypo_paramset', hypo_paramset

    if args.run_mcmc:
        fitter = gf_utils.setup_fitter(args, asimov_paramset)

        ln_prob_eval = partial(
            ln_prob,
            hypo_paramset  = hypo_paramset,
            args           = args,
            fitter         = fitter
        )

        if args.mcmc_seed_type == MCMCSeedType.UNIFORM:
            p0 = mcmc_utils.flat_seed(
                hypo_paramset, nwalkers=args.nwalkers
            )
        elif args.mcmc_seed_type == MCMCSeedType.GAUSSIAN:
            p0 = mcmc_utils.gaussian_seed(
                hypo_paramset, nwalkers=args.nwalkers
            )

        samples = mcmc_utils.mcmc(
            p0       = p0,
            ln_prob  = ln_prob_eval,
            ndim     = n_params,
            nwalkers = args.nwalkers,
            burnin   = args.burnin,
            nsteps   = args.nsteps,
            args     = args,
            threads  = 1
            # TODO(shivesh): broken because you cannot pickle a GolemFitPy object
            # threads      = misc_utils.thread_factors(args.threads)[0]
        )
        mcmc_utils.save_chains(samples, outfile)

    of = outfile[:5]+outfile[5:].replace('data', 'plots')+'_posterior'
    plot_utils.chainer_plot(
        infile       = outfile+'.npy',
        outfile      = of,
        outformat    = ['png'],
        args         = args,
        llh_paramset = hypo_paramset,
        fig_text     = gen_figtext(args, hypo_paramset)
    )

    print "DONE!"


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
