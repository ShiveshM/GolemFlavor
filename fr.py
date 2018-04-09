#! /usr/bin/env python
# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : March 17, 2018

"""
HESE BSM flavour ratio analysis script
"""

from __future__ import absolute_import, division

import argparse
from functools import partial

import numpy as np

import GolemFitPy as gf

from utils import gf as gf_utils
from utils import likelihood as llh_utils
from utils import mcmc as mcmc_utils
from utils import misc as misc_utils
from utils.enums import EnergyDependance, Likelihood, MCMCSeedType, ParamTag
from utils.fr import MASS_EIGENVALUES, normalise_fr
from utils.misc import enum_parse, Param, ParamSet
from utils.plot import plot_argparse, chainer_plot


def define_nuisance():
    """Define the nuisance parameters to scan over with default values,
    ranges and sigma.
    """
    tag = ParamTag.NUISANCE
    nuisance = ParamSet(
        Param(name='convNorm',        value=1., seed=[0.5, 2. ], ranges=[0. , 50.], std=0.3,  tag=tag),
        Param(name='promptNorm',      value=0., seed=[0. , 6. ], ranges=[0. , 50.], std=0.05, tag=tag),
        Param(name='muonNorm',        value=1., seed=[0. , 2. ], ranges=[0. , 50.], std=0.1,  tag=tag),
        Param(name='astroNorm',       value=1., seed=[4. , 10.], ranges=[0. , 50.], std=0.1,  tag=tag),
        Param(name='astroDeltaGamma', value=2., seed=[2. , 3. ], ranges=[-5., 5. ], std=0.1,  tag=tag)
    )
    return nuisance


def nuisance_argparse(parser):
    nuisance_paramset = define_nuisance()
    for parm in nuisance_paramset:
        parser.add_argument(
            '--'+parm.name, type=float, default=parm.value,
            help=parm.name+' to inject'
        )


def get_paramsets(args):
    """Make the paramsets for generating the Asmimov MC sample and also running
    the MCMC.
    """
    asimov_paramset = []
    mcmc_paramset = []
    if args.likelihood == Likelihood.GOLEMFIT:
        nuisance_paramset = define_nuisance()
        asimov_paramset.extend(nuisance_paramset.params)
        mcmc_paramset.extend(nuisance_paramset.params)
        for parm in nuisance_paramset:
            parm.value = args.__getattribute__(parm.name)
    tag = ParamTag.BESTFIT
    asimov_paramset.extend([
        Param(name='astroENorm'  , value=args.measured_ratio[0], ranges=[0., 1.], std=0.2, tag=tag),
        Param(name='astroMuNorm' , value=args.measured_ratio[1], ranges=[0., 1.], std=0.2, tag=tag),
        Param(name='astroTauNorm', value=args.measured_ratio[2], ranges=[0., 1.], std=0.2, tag=tag)
    ])
    asimov_paramset = ParamSet(asimov_paramset)

    if not args.fix_mixing:
        tag = ParamTag.MMANGLES
        mcmc_paramset.extend([
            Param(name='s_12^2', value=0.5, ranges=[0., 1.], std=0.2, tex=r'\tilde{s}_{12}^2', tag=tag),
            Param(name='c_13^4', value=0.5, ranges=[0., 1.], std=0.2, tex=r'\tilde{c}_{13}^4', tag=tag),
            Param(name='s_23^2', value=0.5, ranges=[0., 1.], std=0.2, tex=r'\tilde{s}_{23}^4', tag=tag),
            Param(name='dcp', value=np.pi, ranges=[0., 2*np.pi], std=0.2, tex=r'\tilde{\delta_{CP}}', tag=tag)
        ])
    if not args.fix_scale:
        logLam, scale_region = np.log10(args.scale), np.log10(args.scale_region)
        lL_range = (logLam-scale_region, logLam+scale_region)
        tag = ParamTag.SCALE
        mcmc_paramset.append(
            Param(name='logLam', value=logLam, ranges=lL_range, std=3, tex=r'{\rm log}_{10}\Lambda', tag=tag)
        )
    if not args.fix_source_ratio:
        tag = ParamTag.SRCANGLES
        mcmc_paramset.extend([
            Param(name='s_phi4', value=0.5, ranges=[0., 1.], std=0.2, tex=r'sin^4(\phi)', tag=tag),
            Param(name='c_2psi', value=0.5, ranges=[0., 1.], std=0.2, tex=r'cos(2\psi)', tag=tag)
        ])
    mcmc_paramset = ParamSet(mcmc_paramset)
    # TODO(shivesh): unify
    return asimov_paramset, mcmc_paramset


def process_args(args):
    """Process the input args."""
    if args.fix_mixing and args.fix_source_ratio:
        raise NotImplementedError('Fixed mixing and sfr not implemented')
    if args.fix_mixing and args.fix_scale:
        raise NotImplementedError('Fixed mixing and scale not implemented')

    args.measured_ratio = normalise_fr(args.measured_ratio)
    if args.fix_source_ratio:
        args.source_ratio = normalise_fr(args.source_ratio)

    if args.energy_dependance is EnergyDependance.SPECTRAL:
        args.binning = np.logspace(
            np.log10(args.binning[0]), np.log10(args.binning[1]), args.binning[2]+1
        )

    if not args.fix_scale:
        if args.energy_dependance is EnergyDependance.MONO:
            args.scale = np.power(
                10, np.round(np.log10(MASS_EIGENVALUES[1]/args.energy)) - \
                np.log10(args.energy**(args.dimension-3))
            )
        elif args.energy_dependance is EnergyDependance.SPECTRAL:
            args.scale = np.power(
                10, np.round(
                    np.log10(MASS_EIGENVALUES[1]/np.power(10, np.average(np.log10(args.binning)))) \
                    - np.log10(np.power(10, np.average(np.log10(args.binning)))**(args.dimension-3))
                )
            )
            """Estimate the scale of when to expect the BSM term to contribute."""


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="BSM flavour ratio analysis",
        formatter_class=misc_utils.SortingHelpFormatter,
    )
    parser.add_argument(
        '--measured-ratio', type=int, nargs=3, default=[1, 1, 1],
        help='Set the central value for the measured flavour ratio at IceCube'
    )
    parser.add_argument(
        '--sigma-ratio', type=float, default=0.01,
        help='Set the 1 sigma for the measured flavour ratio for a gaussian LLH'
    )
    parser.add_argument(
        '--fix-source-ratio', type=misc_utils.parse_bool, default='False',
        help='Fix the source flavour ratio'
    )
    parser.add_argument(
        '--energy-dependance', default='spectral', type=partial(enum_parse, c=EnergyDependance),
        choices=EnergyDependance,
        help='Type of energy dependance to use in the BSM fit'
    )
    parser.add_argument(
        '--spectral-index', default=-2, type=int,
        help='Spectral index for spectral energy dependance'
    )
    parser.add_argument(
        '--binning', default=[1e4, 1e7, 10], type=int, nargs=3,
        help='Binning for spectral energy dependance'
    )
    parser.add_argument(
        '--source-ratio', type=int, nargs=3, default=[2, 1, 0],
        help='Set the source flavour ratio for the case when you want to fix it'
    )
    parser.add_argument(
        '--no-bsm', type=misc_utils.parse_bool, default='False',
        help='Turn off BSM terms'
    )
    parser.add_argument(
        '--fix-mixing', type=misc_utils.parse_bool, default='False',
        help='Fix all mixing parameters except one'
    )
    parser.add_argument(
        '--fix-scale', type=misc_utils.parse_bool, default='False',
        help='Fix the new physics scale'
    )
    parser.add_argument(
        '--scale', type=float, default=1e-30,
        help='Set the new physics scale'
    )
    parser.add_argument(
        '--scale-region', type=float, default=1e10,
        help='Set the size of the box to scan for the scale'
    )
    parser.add_argument(
        '--dimension', type=int, default=3,
        help='Set the new physics dimension to consider'
    )
    parser.add_argument(
        '--energy', type=float, default=1000,
        help='Set the energy scale'
    )
    parser.add_argument(
        '--seed', type=int, default=99,
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
    llh_utils.likelihood_argparse(parser)
    mcmc_utils.mcmc_argparse(parser)
    nuisance_argparse(parser)
    gf_utils.gf_argparse(parser)
    plot_argparse(parser)
    return parser.parse_args()


def main():
    args = parse_args()
    process_args(args)
    misc_utils.print_args(args)

    np.random.seed(args.seed)

    asimov_paramset, mcmc_paramset = get_paramsets(args)
    outfile = misc_utils.gen_outfile_name(args)
    print '== {0:<25} = {1}'.format('outfile', outfile)

    if args.run_mcmc:
        if args.likelihood is Likelihood.GOLEMFIT:
            datapaths = gf.DataPaths()
            sparams = gf_utils.steering_params(args)
            npp = gf.NewPhysicsParams()

            fitter = gf.GolemFit(datapaths, sparams, npp)
            gf_utils.set_up_as(fitter, asimov_paramset)
            # fitter.WriteCompact()
        else:
            fitter = None

        ln_prob = partial(
            llh_utils.ln_prob, args=args, fitter=fitter,
            asimov_paramset=asimov_paramset, mcmc_paramset=mcmc_paramset
        )

        ndim = len(mcmc_paramset)
        if args.mcmc_seed_type == MCMCSeedType.UNIFORM:
            p0 = mcmc_utils.flat_seed(
                mcmc_paramset, nwalkers=args.nwalkers
            )
        elif args.mcmc_seed_type == MCMCSeedType.GAUSSIAN:
            p0 = mcmc_utils.gaussian_seed(
                mcmc_paramset, nwalkers=args.nwalkers
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

    print "Making triangle plots"
    chainer_plot(
        infile        = outfile+'.npy',
        outfile       = outfile[:5]+outfile[5:].replace('data', 'plots'),
        outformat     = ['pdf'],
        args          = args,
        mcmc_paramset = mcmc_paramset
    )
    print "DONE!"


main.__doc__ = __doc__


if __name__ == '__main__':
    main()

