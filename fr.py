#! /usr/bin/env python
# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : March 17, 2018

"""
HESE BSM flavour ratio analysis script
"""

from __future__ import absolute_import, division

import os
import argparse
from functools import partial

import numpy as np

import GolemFitPy as gf

from utils import gf as gf_utils
from utils import likelihood as llh_utils
from utils import mcmc as mcmc_utils
from utils import misc as misc_utils
from utils import plot as plot_utils
from utils.enums import EnergyDependance, Likelihood, MCMCSeedType, ParamTag
from utils.fr import MASS_EIGENVALUES, normalise_fr, fr_to_angles
from utils.misc import enum_parse, Param, ParamSet


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
    flavour_angles = fr_to_angles(args.measured_ratio)
    asimov_paramset.extend([
        Param(name='astroFlavorAngle1', value=flavour_angles[0], ranges=[0., 1.], std=0.2, tag=tag),
        Param(name='astroFlavorAngle2', value=flavour_angles[1], ranges=[-1., 1.], std=0.2, tag=tag),
    ])
    asimov_paramset = ParamSet(asimov_paramset)

    if not args.fix_mixing and not args.fix_mixing_almost:
        tag = ParamTag.MMANGLES
        e = 1e-10
        mcmc_paramset.extend([
            Param(name='s_12^2', value=0.5, ranges=[0.+e, 1.-e], std=0.2, tex=r'\tilde{s}_{12}^2', tag=tag),
            Param(name='c_13^4', value=0.5, ranges=[0.+e, 1.-e], std=0.2, tex=r'\tilde{c}_{13}^4', tag=tag),
            Param(name='s_23^2', value=0.5, ranges=[0.+e, 1.-e], std=0.2, tex=r'\tilde{s}_{23}^4', tag=tag),
            Param(name='dcp', value=np.pi, ranges=[0.+e, 2*np.pi-e], std=0.2, tex=r'\tilde{\delta_{CP}}', tag=tag)
        ])
    if args.fix_mixing_almost:
        tag = ParamTag.MMANGLES
        mcmc_paramset.extend([
            Param(name='s_23^2', value=0.5, ranges=[0., 1.], std=0.2, tex=r'\tilde{s}_{23}^4', tag=tag)
        ])
    if not args.fix_scale:
        tag = ParamTag.SCALE
        mcmc_paramset.append(
            Param(name='logLam', value=np.log10(args.scale), ranges=np.log10(args.scale_region), std=3,
                  tex=r'{\rm log}_{10}\Lambda'+plot_utils.get_units(args.dimension), tag=tag)
        )
    if not args.fix_source_ratio:
        tag = ParamTag.SRCANGLES
        mcmc_paramset.extend([
            Param(name='s_phi4', value=0.5, ranges=[0., 1.], std=0.2, tex=r'sin^4(\phi)', tag=tag),
            Param(name='c_2psi', value=0.5, ranges=[0., 1.], std=0.2, tex=r'cos(2\psi)', tag=tag)
        ])
    mcmc_paramset = ParamSet(mcmc_paramset)
    return asimov_paramset, mcmc_paramset


def process_args(args):
    """Process the input args."""
    if args.fix_mixing and args.fix_scale:
        raise NotImplementedError('Fixed mixing and scale not implemented')
    if args.fix_mixing and args.fix_mixing_almost:
        raise NotImplementedError(
            '--fix-mixing and --fix-mixing-almost cannot be used together'
        )
    if args.run_bayes_factor and args.fix_scale:
        raise NotImplementedError(
            '--run-bayes-factor and --fix-scale cannot be used together'
        )

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
                np.log10(args.energy**(args.dimension-3)) + 6
            )
        elif args.energy_dependance is EnergyDependance.SPECTRAL:
            args.scale = np.power(
                10, np.round(
                    np.log10(MASS_EIGENVALUES[1]/np.power(10, np.average(np.log10(args.binning)))) \
                    - np.log10(np.power(10, np.average(np.log10(args.binning)))**(args.dimension-3))
                ) + 6
            )
            """Estimate the scale of when to expect the BSM term to contribute."""
        args.scale_region = (args.scale/args.scale_region, args.scale*args.scale_region)

    if args.bayes_eval_bin.lower() == 'all':
        args.bayes_eval_bin = None
    else:
        args.bayes_eval_bin = int(args.bayes_eval_bin)


def parse_args(args=None):
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
        '--binning', default=[1e4, 1e7, 5], type=float, nargs=3,
        help='Binning for spectral energy dependance'
    )
    parser.add_argument(
        '--run-bayes-factor', type=misc_utils.parse_bool, default='False',
        help='Make the bayes factor plot for the scale'
    )
    parser.add_argument(
        '--run-angles-limit', type=misc_utils.parse_bool, default='False',
        help='Make the limit vs BSM angles plot'
    )
    parser.add_argument(
        '--run-angles-correlation', type=misc_utils.parse_bool, default='False',
        help='Make the limit vs BSM angles plot'
    )
    parser.add_argument(
        '--bayes-bins', type=int, default=10,
        help='Number of bins for the Bayes factor plot'
    )
    parser.add_argument(
        '--bayes-eval-bin', type=str, default='all',
        help='Which bin to evalaute for Bayes factor plot'
    )
    parser.add_argument(
        '--bayes-live-points', type=int, default=400,
        help='Number of live points for MultiNest evaluations'
    )
    parser.add_argument(
        '--bayes-tolerance', type=float, default=0.5,
        help='Tolerance for MultiNest'
    )
    parser.add_argument(
        '--bayes-output', type=str, default='./mnrun/',
        help='Folder to store MultiNest evaluations'
    )
    parser.add_argument(
        '--angles-lim-output', type=str, default='./mnrun/',
        help='Folder to store MultiNest evaluations'
    )
    parser.add_argument(
        '--angles-corr-output', type=str, default='./mnrun/',
        help='Folder to store MultiNest evaluations'
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
        help='Fix all mixing parameters to bi-maximal mixing'
    )
    parser.add_argument(
        '--fix-mixing-almost', type=misc_utils.parse_bool, default='False',
        help='Fix all mixing parameters except s23'
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
        '--scale-region', type=float, default=1e7,
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
    gf_utils.gf_argparse(parser)
    plot_utils.plot_argparse(parser)
    nuisance_argparse(parser)
    if args is None: return parser.parse_args()
    else: return parser.parse_args(args.split())


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
            fitter = gf_utils.setup_fitter(args, asimov_paramset)
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

    plot_utils.chainer_plot(
        infile        = outfile+'.npy',
        outfile       = outfile[:5]+outfile[5:].replace('data', 'plots'),
        outformat     = ['pdf'],
        args          = args,
        mcmc_paramset = mcmc_paramset
    )

    out = args.bayes_output+'/fr_evidence' + misc_utils.gen_identifier(args)
    scan_scales = np.linspace(
        np.log10(args.scale_region[0]), np.log10(args.scale_region[1]), args.bayes_bins
    )
    if args.run_bayes_factor:
        import pymultinest

        if not args.run_mcmc and args.likelihood is Likelihood.GOLEMFIT:
            fitter = gf_utils.setup_fitter(args, asimov_paramset)
        else: fitter = None

        p = mcmc_paramset.from_tag(ParamTag.SCALE, invert=True)
        n_params = len(p)
        prior_ranges = p.seeds

        def CubePrior(cube, ndim, nparams):
            # default are uniform priors
            return ;

        arr = []
        for s_idx, sc in enumerate(scan_scales):
            if args.bayes_eval_bin is not None:
                if args.bayes_eval_bin == s_idx:
                    out += '_scale_{0:.0E}'.format(np.power(10, sc))
                else: continue

            print '== SCALE = {0:.0E}'.format(np.power(10, sc))
            theta = np.zeros(n_params)
            def lnProb(cube, ndim, nparams):
                for i in range(ndim):
                    prange = prior_ranges[i][1] - prior_ranges[i][0]
                    theta[i] = prange*cube[i] + prior_ranges[i][0]
                theta_ = np.array(theta.tolist() + [sc])
                # print 'mcmc_paramset', mcmc_paramset
                return llh_utils.triangle_llh(
                    theta=theta_,
                    args=args,
                    asimov_paramset=asimov_paramset,
                    mcmc_paramset=mcmc_paramset,
                    fitter=fitter
                )

            prefix = 'mnrun/' + os.path.basename(out) + '_'
            if args.bayes_eval_bin is not None:
                prefix += '{0}_'.format(s_idx)
            print 'begin running evidence calculation for {0}'.format(prefix)
            result = pymultinest.run(
                LogLikelihood=lnProb,
                Prior=CubePrior,
                n_dims=n_params,
                importance_nested_sampling=True,
                n_live_points=args.bayes_live_points,
                evidence_tolerance=args.bayes_tolerance,
                outputfiles_basename=prefix,
                resume=False,
                verbose=True
            )

            analyzer = pymultinest.Analyzer(
                outputfiles_basename=prefix, n_params=n_params
            )
            a_lnZ = analyzer.get_stats()['global evidence']
            print 'Evidence = {0}'.format(a_lnZ)
            arr.append([sc, a_lnZ])

        misc_utils.make_dir(out)
        np.save(out+'.npy', np.array(arr))

    dirname = os.path.dirname(out)
    plot_utils.bayes_factor_plot(
        dirname=dirname, outfile=out, outformat=['png'], args=args
    )

    out = args.angles_lim_output+'/fr_an_evidence' + misc_utils.gen_identifier(args)
    if args.run_angles_limit:
        import pymultinest

        scenarios = [
            [np.sin(np.pi/2.)**2, 0, 0, 0],
            [0, np.cos(np.pi/2.)**4, 0, 0],
            [0, 0, np.sin(np.pi/2.)**2, 0],
        ]
        p = mcmc_paramset.from_tag([ParamTag.SCALE, ParamTag.MMANGLES], invert=True)
        n_params = len(p)
        prior_ranges = p.seeds

        if not args.run_mcmc and args.likelihood is Likelihood.GOLEMFIT:
            fitter = gf_utils.setup_fitter(args, asimov_paramset)
        else: fitter = None

        def CubePrior(cube, ndim, nparams):
            # default are uniform priors
            return ;

        if args.bayes_eval_bin is not None:
            data = np.zeros((len(scenarios), 1, 2))
        else: data = np.zeros((len(scenarios), args.bayes_bins, 2))
        mm_angles = mcmc_paramset.from_tag(ParamTag.MMANGLES)
        sc_angles = mcmc_paramset.from_tag(ParamTag.SCALE)[0]
        for idx, scen in enumerate(scenarios):
            scales, evidences = [], []
            for yidx, an in enumerate(mm_angles):
                an.value = scen[yidx]
            for s_idx, sc in enumerate(scan_scales):
                if args.bayes_eval_bin is not None:
                    if args.bayes_eval_bin == s_idx:
                        if idx == 0:
                            out += '_scale_{0:.0E}'.format(np.power(10, sc))
                    else: continue

                print '== SCALE = {0:.0E}'.format(np.power(10, sc))
                sc_angles.value = sc
                def lnProb(cube, ndim, nparams):
                    for i in range(ndim):
                        prange = prior_ranges[i][1] - prior_ranges[i][0]
                        p[i].value = prange*cube[i] + prior_ranges[i][0]
                    for name in p.names:
                        mcmc_paramset[name].value = p[name].value
                    theta = mcmc_paramset.values
                    # print 'theta', theta
                    # print 'mcmc_paramset', mcmc_paramset
                    return llh_utils.triangle_llh(
                        theta=theta,
                        args=args,
                        asimov_paramset=asimov_paramset,
                        mcmc_paramset=mcmc_paramset,
                        fitter=fitter
                    )
                prefix = 'mnrun/' + os.path.basename(out) + '_'
                if args.bayes_eval_bin is not None:
                    prefix += '{0}_{1}_'.format(idx, s_idx)
                print 'begin running evidence calculation for {0}'.format(prefix)
                result = pymultinest.run(
                    LogLikelihood=lnProb,
                    Prior=CubePrior,
                    n_dims=n_params,
                    importance_nested_sampling=True,
                    n_live_points=args.bayes_live_points,
                    evidence_tolerance=args.bayes_tolerance,
                    outputfiles_basename=prefix,
                    resume=False,
                    verbose=True
                )

                analyzer = pymultinest.Analyzer(outputfiles_basename=prefix, n_params=n_params)
                a_lnZ = analyzer.get_stats()['global evidence']
                print 'Evidence = {0}'.format(a_lnZ)
                scales.append(sc)
                evidences.append(a_lnZ)

            for i, d in enumerate(evidences):
                data[idx][i][0] = scales[i]
                data[idx][i][1] = d

        misc_utils.make_dir(out)
        print 'saving to {0}.npy'.format(out)
        np.save(out+'.npy', np.array(data))

    dirname = os.path.dirname(out)
    plot_utils.plot_BSM_angles_limit(
        dirname=dirname, outfile=outfile, outformat=['png'],
        args=args, bayesian=True
    )

    out = args.angles_corr_output+'/fr_co_evidence' + misc_utils.gen_identifier(args)
    if args.run_angles_correlation:
        if args.bayes_eval_bin is None: assert 0
        import pymultinest

        scenarios = [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
        ]
        nuisance = mcmc_paramset.from_tag(ParamTag.NUISANCE)
        mm_angles = mcmc_paramset.from_tag(ParamTag.MMANGLES)
        sc_angles = mcmc_paramset.from_tag(ParamTag.SCALE)[0]

        if not args.run_mcmc and args.likelihood is Likelihood.GOLEMFIT:
            fitter = gf_utils.setup_fitter(args, asimov_paramset)
        else: fitter = None

        def CubePrior(cube, ndim, nparams):
            # default are uniform priors
            return ;

        scan_angles = np.linspace(0, 1, args.bayes_bins)

        if args.bayes_eval_bin is not None:
            data = np.zeros((len(scenarios), 1, 1, 3))
        else: data = np.zeros((len(scenarios), scale_bins, angle_bins, 3))
        for idx, scen in enumerate(scenarios):
            for an in mm_angles:
                an.value = 0
            keep = mcmc_paramset.from_tag(ParamTag.MMANGLES)[idx]
            q = ParamSet(nuisance.params + [x for x in mm_angles if x.name != keep.name])
            n_params = len(q)
            prior_ranges = q.seeds

            scales, angles, evidences = [], [], []
            for s_idx, sc in enumerate(scan_scales):
                for a_idx, an in enumerate(scan_angles):
                    index = s_idx*args.bayes_bins + a_idx
                    if args.bayes_eval_bin is not None:
                        if args.bayes_eval_bin == index:
                            if idx == 0:
                                out += '_scale_{0:.0E}'.format(np.power(10, sc))
                        else: continue

                    print '== SCALE = {0:.0E}'.format(np.power(10, sc))
                    print '== ANGLE = {0:.2f}'.format(an)
                    sc_angles.value = sc
                    keep.value = an
                    def lnProb(cube, ndim, nparams):
                        for i in range(ndim-1):
                            prange = prior_ranges[i][1] - prior_ranges[i][0]
                            q[i].value = prange*cube[i] + prior_ranges[i][0]
                        for name in q.names:
                            mcmc_paramset[name].value = q[name].value
                        theta = mcmc_paramset.values
                        # print 'theta', theta
                        # print 'mcmc_paramset', mcmc_paramset
                        return llh_utils.triangle_llh(
                            theta=theta,
                            args=args,
                            asimov_paramset=asimov_paramset,
                            mcmc_paramset=mcmc_paramset,
                            fitter=fitter
                        )
                    prefix = 'mnrun/' + os.path.basename(out) + '_'
                    if args.bayes_eval_bin is not None:
                        prefix += '{0}_{1}_{2}'.format(idx, s_idx, a_idx)

                    print 'begin running evidence calculation for {0}'.format(prefix)
                    result = pymultinest.run(
                        LogLikelihood=lnProb,
                        Prior=CubePrior,
                        n_dims=n_params,
                        importance_nested_sampling=True,
                        n_live_points=args.bayes_live_points,
                        evidence_tolerance=args.bayes_tolerance,
                        outputfiles_basename=prefix,
                        resume=False,
                        verbose=True
                    )

                    analyzer = pymultinest.Analyzer(outputfiles_basename=prefix, n_params=n_params)
                    a_lnZ = analyzer.get_stats()['global evidence']
                    print 'Evidence = {0}'.format(a_lnZ)
                    scales.append(sc)
                    angles.append(an)
                    evidences.append(a_lnZ)

                for i, d in enumerate(evidences):
                    data[idx][i][i][0] = scales[i]
                    data[idx][i][i][1] = angles[i]
                    data[idx][i][i][2] = d

        misc_utils.make_dir(out)
        print 'saving to {0}.npy'.format(out)
        np.save(out+'.npy', np.array(data))

    print "DONE!"


main.__doc__ = __doc__


if __name__ == '__main__':
    main()

