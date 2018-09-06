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
from utils.enums import EnergyDependance, Likelihood, MixingScenario
from utils.enums import MCMCSeedType, ParamTag, PriorsCateg
from utils.param import Param, ParamSet, get_paramsets


def define_nuisance():
    """Define the nuisance parameters."""
    tag = ParamTag.SM_ANGLES
    nuisance = []
    g_prior = PriorsCateg.GAUSSIAN
    hg_prior = PriorsCateg.HALFGAUSS
    e = 1e-9
    nuisance.extend([
        Param(name='s_12_2', value=0.307,            seed=[0.26, 0.35],     ranges=[0., 1.],      std=0.013,   tex=r's_{12}^2', prior=g_prior,  tag=tag),
        Param(name='c_13_4', value=(1-(0.02206))**2, seed=[0.950, 0.961],   ranges=[0., 1.],      std=0.00147, tex=r'c_{13}^4', prior=g_prior, tag=tag),
        Param(name='s_23_2', value=0.538,            seed=[0.31, 0.75],     ranges=[0., 1.],      std=0.069,   tex=r's_{23}^2', prior=g_prior,  tag=tag),
        Param(name='dcp',    value=4.08404,          seed=[0+e, 2*np.pi-e], ranges=[0., 2*np.pi], std=2.0,     tex=r'\delta_{CP}', tag=tag),
        Param(
            name='m21_2', value=7.40E-23, seed=[7.2E-23, 7.6E-23], ranges=[6.80E-23, 8.02E-23],
            std=2.1E-24, tex=r'\Delta m_{21}^2{\rm GeV}^{-2}', prior=g_prior, tag=tag
        ),
        Param(
            name='m3x_2', value=2.494E-21, seed=[2.46E-21, 2.53E-21], ranges=[2.399E-21, 2.593E-21],
            std=3.3E-23, tex=r'\Delta m_{3x}^2{\rm GeV}^{-2}', prior=g_prior, tag=tag
        )
    ])
    tag = ParamTag.NUISANCE
    nuisance.extend([
        Param(name='convNorm',        value=1.,  seed=[0.5, 2. ], ranges=[0. , 50.], std=0.3,  tag=tag),
        Param(name='promptNorm',      value=0.,  seed=[0. , 6. ], ranges=[0. , 50.], std=0.05, tag=tag),
        Param(name='muonNorm',        value=1.,  seed=[0.1, 2. ], ranges=[0. , 50.], std=0.1,  tag=tag),
        Param(name='astroNorm',       value=6.9, seed=[0.1, 10.], ranges=[0. , 50.], std=0.1,  tag=tag),
        Param(name='astroDeltaGamma', value=2.5, seed=[1. , 3. ], ranges=[-5., 5. ], std=0.1,  tag=tag)
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
    if args.fix_mixing is not MixingScenario.NONE and args.fix_scale:
        raise NotImplementedError('Fixed mixing and scale not implemented')
    if args.fix_mixing is not MixingScenario.NONE and args.fix_mixing_almost:
        raise NotImplementedError(
            '--fix-mixing and --fix-mixing-almost cannot be used together'
        )

    args.measured_ratio = fr_utils.normalise_fr(args.measured_ratio)
    if args.fix_source_ratio:
        args.source_ratio = fr_utils.normalise_fr(args.source_ratio)

    if args.energy_dependance is EnergyDependance.SPECTRAL:
        args.binning = np.logspace(
            np.log10(args.binning[0]), np.log10(args.binning[1]), args.binning[2]+1
        )

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

    asimov_paramset, llh_paramset = get_paramsets(args, define_nuisance())
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
