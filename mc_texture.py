#! /usr/bin/env python
# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : April 25, 2019

"""
Sample points for a specific scenario
"""

from __future__ import absolute_import, division

import argparse
from copy import deepcopy
from functools import partial

import numpy as np

from utils import fr as fr_utils
from utils import llh as llh_utils
from utils import mcmc as mcmc_utils
from utils import misc as misc_utils
from utils import plot as plot_utils
from utils.enums import MCMCSeedType, ParamTag, PriorsCateg, Texture
from utils.param import Param, ParamSet


def define_nuisance():
    """Define the nuisance parameters."""
    tag = ParamTag.SM_ANGLES
    nuisance = []
    g_prior = PriorsCateg.GAUSSIAN
    lg_prior = PriorsCateg.LIMITEDGAUSS
    e = 1e-9
    nuisance.extend([
        Param(name='s_12_2', value=0.307,            seed=[0.26, 0.35],     ranges=[0., 1.],      std=0.013,   tex=r's_{12}^2', prior=lg_prior,  tag=tag),
        Param(name='c_13_4', value=(1-(0.02206))**2, seed=[0.950, 0.961],   ranges=[0., 1.],      std=0.00147, tex=r'c_{13}^4', prior=lg_prior,  tag=tag),
        Param(name='s_23_2', value=0.538,            seed=[0.31, 0.75],     ranges=[0., 1.],      std=0.069,   tex=r's_{23}^2', prior=lg_prior,  tag=tag),
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
    return ParamSet(nuisance)


def get_paramsets(args, nuisance_paramset):
    """Make the paramsets for generating the Asmimov MC sample and also running
    the MCMC.
    """
    asimov_paramset = []
    llh_paramset = []

    llh_paramset.extend(
        [x for x in nuisance_paramset.from_tag(ParamTag.SM_ANGLES)]
    )

    for parm in llh_paramset:
        parm.value = args.__getattribute__(parm.name)

    boundaries = fr_utils.SCALE_BOUNDARIES[args.dimension]
    tag = ParamTag.SCALE
    llh_paramset.append(
        Param(
            name='logLam', value=np.mean(boundaries), ranges=boundaries, std=3,
            tex=r'{\rm log}_{10}\left (\Lambda^{-1}' + \
                misc_utils.get_units(args.dimension)+r'\right )',
            tag=tag
        )
    )
    llh_paramset = ParamSet(llh_paramset)

    tag = ParamTag.BESTFIT
    flavour_angles = fr_utils.fr_to_angles([1, 1, 1])
    asimov_paramset.extend([
        Param(name='astroFlavorAngle1', value=flavour_angles[0], ranges=[ 0., 1.], std=0.2, tag=tag),
        Param(name='astroFlavorAngle2', value=flavour_angles[1], ranges=[-1., 1.], std=0.2, tag=tag),
    ])
    asimov_paramset = ParamSet(asimov_paramset)

    return asimov_paramset, llh_paramset


def nuisance_argparse(parser):
    nuisance = define_nuisance()
    for parm in nuisance:
        parser.add_argument(
            '--'+parm.name, type=float, default=parm.value,
            help=parm.name+' to inject'
        )


def process_args(args):
    """Process the input args."""
    if args.texture is Texture.NONE:
        raise ValueError('Must assume a BSM texture')
    args.source_ratio = fr_utils.normalise_fr(args.source_ratio)

    args.binning = np.logspace(
        np.log10(args.binning[0]), np.log10(args.binning[1]), args.binning[2]+1
    )


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
        '--spectral-index', type=float, default='-2',
        help='Astro spectral index'
    )
    parser.add_argument(
        '--datadir', type=str, default='./untitled',
        help='Path to store chains'
    )
    fr_utils.fr_argparse(parser)
    mcmc_utils.mcmc_argparse(parser)
    nuisance_argparse(parser)
    misc_utils.remove_option(parser, 'injected_ratio')
    if args is None: return parser.parse_args()
    else: return parser.parse_args(args.split())


def gen_identifier(args):
    f = '_DIM{0}'.format(args.dimension)
    f += '_sfr_' + misc_utils.solve_ratio(args.source_ratio)
    f += '_{0}'.format(misc_utils.str_enum(args.texture))
    return f


def gen_figtext(args):
    """Generate the figure text."""
    t = r'$'
    t += r'{\rm Source\:flavour\:ratio}'+r'\:=\:({0})'.format(
        misc_utils.solve_ratio(args.source_ratio).replace('_', ':')
    )
    t += '$\n' + r'${\rm Texture}'+r' = {0}'.format(
        misc_utils.str_enum(args.texture)
    )
    t += '$\n' + r'${\rm Dimension}'+r' = {0}$'.format(args.dimension)
    return t


def triangle_llh(theta, args, llh_paramset):
    """Log likelihood function for a given theta."""
    if len(theta) != len(llh_paramset):
        raise AssertionError(
            'Dimensions of scan is not the same as the input '
            'params\ntheta={0}\nparamset]{1}'.format(theta, llh_paramset)
        )
    for idx, param in enumerate(llh_paramset):
        param.value = theta[idx]

    return 1. # Flat LLH


def ln_prob(theta, args, llh_paramset):
    dc_llh_paramset = deepcopy(llh_paramset)
    lp = llh_utils.lnprior(theta, paramset=dc_llh_paramset)
    if not np.isfinite(lp):
        return -np.inf
    return lp + triangle_llh(
        theta,
        args          = args,
        llh_paramset  = dc_llh_paramset,
    )


def main():
    args = parse_args()
    process_args(args)
    misc_utils.print_args(args)

    if args.seed is not None:
        np.random.seed(args.seed)

    asimov_paramset, llh_paramset = get_paramsets(args, define_nuisance())

    prefix = ''
    outfile = args.datadir + '/mc_texture' + prefix + gen_identifier(args)
    print '== {0:<25} = {1}'.format('outfile', outfile)

    print 'asimov_paramset', asimov_paramset
    print 'llh_paramset', llh_paramset

    if args.run_mcmc:
        ln_prob_eval = partial(
            ln_prob,
            llh_paramset  = llh_paramset,
            args           = args,
        )

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
            ln_prob  = ln_prob_eval,
            ndim     = len(llh_paramset),
            nwalkers = args.nwalkers,
            burnin   = args.burnin,
            nsteps   = args.nsteps,
            args     = args,
            threads  = args.threads
        )

        frs = np.array(
            map(lambda x: fr_utils.flux_averaged_BSMu(
                x, args, args.spectral_index, llh_paramset
            ), samples),
            dtype=float
        )
        mcmc_utils.save_chains(frs, outfile)

    of = outfile[:5]+outfile[5:].replace('data', 'plots')+'_posterior'
    plot_utils.chainer_plot(
        infile       = outfile+'.npy',
        outfile      = of,
        outformat    = ['png'],
        args         = args,
        llh_paramset = llh_paramset,
        fig_text     = gen_figtext(args, llh_paramset)
    )
    print "DONE!"


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
