#! /usr/bin/env python
# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : March 17, 2018

"""
emcee sample SM points for (x, 1-x, 0)
"""

from __future__ import absolute_import, division

import argparse
from copy import deepcopy
from functools import partial

import numpy as np

from golemflavor import fr as fr_utils
from golemflavor import llh as llh_utils
from golemflavor import mcmc as mcmc_utils
from golemflavor import misc as misc_utils
from golemflavor import plot as plot_utils
from golemflavor.enums import MCMCSeedType, ParamTag, PriorsCateg
from golemflavor.param import Param, ParamSet


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
    ])
    tag = ParamTag.SRCANGLES
    nuisance.extend([
        Param(name='astroX', value=0.5, seed=[0., 1.], ranges=[0., 1.], std=0.1, tex=r'x', tag=tag)
    ])
    return ParamSet(nuisance)


def get_paramsets(args, nuisance_paramset):
    """Make the paramsets for generating the Asmimov MC sample and also running
    the MCMC.
    """
    asimov_paramset = []
    hypo_paramset = []

    hypo_paramset.extend(
        [x for x in nuisance_paramset.from_tag((
            ParamTag.SM_ANGLES, ParamTag.SRCANGLES
        ))]
    )

    for parm in hypo_paramset:
        parm.value = args.__getattribute__(parm.name)

    hypo_paramset = ParamSet(hypo_paramset)

    tag = ParamTag.BESTFIT
    asimov_paramset.extend([
        Param(name='astroFlavorAngle1', value=0., ranges=[ 0., 1.], std=0.2, tag=tag),
        Param(name='astroFlavorAngle2', value=0., ranges=[-1., 1.], std=0.2, tag=tag),
    ])
    asimov_paramset = ParamSet(asimov_paramset)

    return asimov_paramset, hypo_paramset


def nuisance_argparse(parser):
    nuisance = define_nuisance()
    for parm in nuisance:
        parser.add_argument(
            '--'+parm.name, type=float, default=parm.value,
            help=parm.name+' to inject'
        )


def process_args(args):
    """Process the input args."""
    pass


def parse_args(args=None):
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="BSM flavor ratio analysis",
        formatter_class=misc_utils.SortingHelpFormatter,
    )
    parser.add_argument(
        '--seed', type=misc_utils.seed_parse, default='26',
        help='Set the random seed value'
    )
    parser.add_argument(
        '--threads', type=misc_utils.thread_type, default='1',
        help='Set the number of threads to use (int or "max")'
    )
    parser.add_argument(
        '--datadir', type=str, default='./untitled',
        help='Path to store chains'
    )
    mcmc_utils.mcmc_argparse(parser)
    nuisance_argparse(parser)
    misc_utils.remove_option(parser, 'injected_ratio')
    misc_utils.remove_option(parser, 'plot_angles')
    misc_utils.remove_option(parser, 'plot_elements')
    if args is None: return parser.parse_args()
    else: return parser.parse_args(args.split())


def triangle_llh(theta, args, hypo_paramset):
    """Log likelihood function for a given theta."""
    if len(theta) != len(hypo_paramset):
        raise AssertionError(
            'Dimensions of scan is not the same as the input '
            'params\ntheta={0}\nparamset]{1}'.format(theta, hypo_paramset)
        )
    for idx, param in enumerate(hypo_paramset):
        param.value = theta[idx]

    return 1. # Flat LLH


def ln_prob(theta, args, hypo_paramset):
    dc_hypo_paramset = deepcopy(hypo_paramset)
    lp = llh_utils.lnprior(theta, paramset=dc_hypo_paramset)
    if not np.isfinite(lp):
        return -np.inf
    return lp + triangle_llh(
        theta,
        args           = args,
        hypo_paramset  = dc_hypo_paramset,
    )


def main():
    args = parse_args()
    process_args(args)
    misc_utils.print_args(args)

    if args.seed is not None:
        np.random.seed(args.seed)

    asimov_paramset, hypo_paramset = get_paramsets(args, define_nuisance())

    prefix = ''
    outfile = args.datadir + '/mc_x' + prefix
    print('== {0:<25} = {1}'.format('outfile', outfile))

    print('asimov_paramset', asimov_paramset)
    print('hypo_paramset', hypo_paramset)

    if args.run_mcmc:
        ln_prob_eval = partial(
            ln_prob,
            hypo_paramset  = hypo_paramset,
            args           = args,
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
            ndim     = len(hypo_paramset),
            nwalkers = args.nwalkers,
            burnin   = args.burnin,
            nsteps   = args.nsteps,
            args     = args,
            threads  = args.threads
        )

        nsamples = len(samples)
        srcs = [fr_utils.normalize_fr((x, 1-x, 0)) for x in samples.T[-1]]
        mmxs = map(fr_utils.angles_to_u, samples.T[:-1].T)
        frs = np.array(
            [fr_utils.u_to_fr(srcs[i], mmxs[i]) for i in range(nsamples)],
            dtype=np.float64
        )
        mcmc_utils.save_chains(frs, outfile)

    print("DONE!")


main.__doc__ = __doc__


if __name__ == '__main__':
    main()

