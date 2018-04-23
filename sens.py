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
import numpy.ma as ma

from utils import fr as fr_utils
from utils import gf as gf_utils
from utils import likelihood as llh_utils
from utils import misc as misc_utils
from utils import plot as plot_utils
from utils.enums import EnergyDependance, Likelihood, ParamTag
from utils.enums import SensitivityCateg, StatCateg
from utils.fr import estimate_scale, normalise_fr
from utils.misc import enum_parse, parse_bool
from utils.param import Param, ParamSet, get_paramsets

from utils import multinest as mn_utils


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
    if args.mn_run and args.fix_scale:
        raise NotImplementedError(
            '--mn-run and --fix-scale cannot be used together'
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

    if args.mn_eval_bin.lower() == 'all':
        args.mn_eval_bin = None
    else:
        args.mn_eval_bin = int(args.mn_eval_bin)


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
        '--run-method', default='full',
        type=partial(enum_parse, c=SensitivityCateg), choices=SensitivityCateg,
        help='Choose which type of sensivity plot to make'
    )
    parser.add_argument(
        '--stat-method', default='bayesian',
        type=partial(enum_parse, c=StatCateg), choices=StatCateg,
        help='Statistical method to employ'
    )
    fr_utils.fr_argparse(parser)
    gf_utils.gf_argparse(parser)
    llh_utils.likelihood_argparse(parser)
    mn_utils.mn_argparse(parser)
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
    scale = llh_paramset.from_tag(ParamTag.SCALE)[0]
    outfile = misc_utils.gen_outfile_name(args)
    print '== {0:<25} = {1}'.format('outfile', outfile)

    mmangles = llh_paramset.from_tag(ParamTag.MMANGLES)
    if args.run_method is SensitivityCateg.FULL:
        mn_paramset_arr = [llh_paramset.from_tag(ParamTag.SCALE, invert=True)]
    elif args.run_method is SensitivityCateg.FIXED_ANGLE or \
            args.run_method is SensitivityCateg.CORR_ANGLE:
        nscale_pmset = llh_paramset.from_tag(ParamTag.SCALE, invert=True)
        mn_paramset_arr = []
        for x in xrange(3):
            mn_paramset_arr.append(
                ParamSet([prms for prms in nscale_pmset
                          if mmangles[x].name != prms.name])
            )

    out = args.outfile+'/{0}/{1}/fr_evidence'.format(args.stat_method, args.run_method) \
        + misc_utils.gen_identifier(args)
    if args.mn_run:
        if args.likelihood is Likelihood.GOLEMFIT:
            fitter = gf_utils.setup_fitter(args, asimov_paramset)
        else: fitter = None

        scan_scales = np.linspace(
            np.log10(args.scale_region[0]), np.log10(args.scale_region[1]), args.mn_bins
        )
        if args.run_method is SensitivityCateg.CORR_ANGLE:
            scan_angles = np.linspace(0, 1, eval_dim)
        else: scan_angles = np.array([0])
        print 'scan_scales', scan_scales
        print 'scan_angles', scan_angles

        if args.mn_eval_bin is None:
            eval_dim = args.mn_bins
        else: eval_dim = 1

        if args.run_method is SensitivityCateg.FULL:
            evidence_arr = np.full((eval_dim, 2), np.nan)
        elif args.run_method is SensitivityCateg.FIXED_ANGLE:
            evidence_arr = np.full((len(mn_paramset_arr), eval_dim, 2), np.nan)
        elif args.run_method is SensitivityCateg.CORR_ANGLE:
            evidence_arr = np.full((len(mn_paramset_arr), eval_dim, eval_dim, 3), np.nan)


        for idx_scen, mn_paramset in enumerate(mn_paramset_arr):
            print '|||| SCENARIO = {0}'.format(idx_scen)
            for x in mmangles: x.value = 0.
            if args.run_method is SensitivityCateg.FIXED_ANGLE:
                if idx_scen == 0 or idx_scen == 2:
                    mmangles[idx_scen].value = np.sin(np.pi/2.)**2
                    """s_12^2 or s_23^2"""
                elif idx_scen == 1:
                    mmangles[idx_scen].value = np.cos(np.pi/2.)**4
                    """c_13^4"""

            for idx_an, an in enumerate(scan_angles):
                if args.run_method is SensitivityCateg.CORR_ANGLE:
                    print '|||| ANGLE = {0:<04.2}'.format(an)
                    mmangles[idx_an].value = an

                for idx_sc, sc in enumerate(scan_scales):
                    if args.mn_eval_bin is not None:
                        if idx_sc == args.mn_eval_bin:
                            out += '_scale_{0:.0E}'.format(np.power(10, sc))
                            if args.run_method is SensitivityCateg.CORR_ANGLE:
                                out += '_angle_{0:>03}'.format(int(an*100))
                            break
                        else: continue

                    print '|||| SCALE = {0:.0E}'.format(np.power(10, sc))
                    scale.value = sc
                    try:
                        a_lnZ = mn_utils.mn_evidence(
                            mn_paramset     = mn_paramset,
                            llh_paramset    = llh_paramset,
                            asimov_paramset = asimov_paramset,
                            args            = args,
                            fitter          = fitter
                        )
                    except:
                        print 'Failed run, continuing'
                        continue
                    print '## Evidence = {0}'.format(a_lnZ)
                    if args.run_method is SensitivityCateg.FULL:
                        evidence_arr[idx_sc] = np.array([sc, a_lnZ])
                    elif args.run_method is SensitivityCateg.FIXED_ANGLE:
                        evidence_arr[idx_scen][idx_sc] = np.array([sc, a_lnZ])
                    elif args.run_method is SensitivityCateg.CORR_ANGLE:
                        evidence_arr[idx_scen][idx_an][idx_sc] = np.array([an, sc, a_lnZ])

        misc_utils.make_dir(out)
        print 'Saving to {0}'.format(out+'.npy')
        np.save(out+'.npy', np.array(evidence_arr))

    if args.plot_multinest:
        if args.mn_run: raw = evidence_arr
        else: raw = np.load(out+'.npy')
        data = ma.masked_invalid(raw, 0)

        basename = os.path.dirname(out) + '/mnrun/' + os.path.basename(out)
        baseoutfile = basename[:5]+basename[5:].replace('data', 'plots')
        if args.run_method is SensitivityCateg.FULL:
            plot_utils.plot_multinest(
                data        = data,
                outfile     = baseoutfile,
                outformat   = ['png'],
                args        = args,
                scale_param = scale
            )
        elif args.run_method is SensitivityCateg.FIXED_ANGLE:
            for idx_scen in xrange(len(mn_paramset_arr)):
                print '|||| SCENARIO = {0}'.format(idx_scen)
                outfile = baseoutfile + '_SCEN{0}'.format(idx_scen)
                if idx_scen == 0: label = r'$\mathcal{O}_{12}=\frac{1}{2}$'
                elif idx_scen == 1: label = r'$\mathcal{O}_{13}=\frac{1}{2}$'
                elif idx_scen == 2: label = r'$\mathcal{O}_{23}=\frac{1}{2}$'
                plot_utils.plot_multinest(
                    data        = data[idx_scen],
                    outfile     = outfile,
                    outformat   = ['png'],
                    args        = args,
                    scale_param = scale,
                    label       = label
                )
        elif args.run_method is SensitivityCateg.CORR_ANGLE:
            for idx_scen in xrange(len(mn_paramset_arr)):
                print '|||| SCENARIO = {0}'.format(idx_scen)
                basescenoutfile = baseoutfile + '_SCEN{0}'.format(idx_scen)
                if idx_scen == 0: label = r'$\mathcal{O}_{12}='
                elif idx_scen == 1: label = r'$\mathcal{O}_{13}='
                elif idx_scen == 2: label = r'$\mathcal{O}_{23}='
                for idx_an, an in enumerate(scan_angles):
                    print '|||| ANGLE = {0:<04.2}'.format(an)
                    outfile = basescenoutfile + '_ANGLE{0}'.format(idx_an)
                    label += r'{0:<04.2}$'.format(an)
                    plot_utils.plot_multinest(
                        data        = data[idx_scen][idx_an][:,1:],
                        outfile     = outfile,
                        outformat   = ['png'],
                        args        = args,
                        scale_param = scale,
                        label       = label
                    )

    # dirname = os.path.dirname(out)
    # plot_utils.bayes_factor_plot(
    #     dirname=dirname, outfile=out, outformat=['png'], args=args
    # )

    # out = args.angles_lim_output+'/fr_an_evidence' + misc_utils.gen_identifier(args)
    # if args.run_angles_limit:
    #     import pymultinest

    #     scenarios = [
    #         [np.sin(np.pi/2.)**2, 0, 0, 0],
    #         [0, np.cos(np.pi/2.)**4, 0, 0],
    #         [0, 0, np.sin(np.pi/2.)**2, 0],
    #     ]
    #     p = llh_paramset.from_tag([ParamTag.SCALE, ParamTag.MMANGLES], invert=True)
    #     n_params = len(p)
    #     prior_ranges = p.seeds

    #     if not args.run_llh and args.likelihood is Likelihood.GOLEMFIT:
    #         fitter = gf_utils.setup_fitter(args, asimov_paramset)
    #     else: fitter = None

    #     def CubePrior(cube, ndim, nparams):
    #         # default are uniform priors
    #         return ;

    #     if args.bayes_eval_bin is not None:
    #         data = np.zeros((len(scenarios), 1, 2))
    #     else: data = np.zeros((len(scenarios), args.bayes_bins, 2))
    #     mm_angles = llh_paramset.from_tag(ParamTag.MMANGLES)
    #     sc_angles = llh_paramset.from_tag(ParamTag.SCALE)[0]
    #     for idx, scen in enumerate(scenarios):
    #         scales, evidences = [], []
    #         for yidx, an in enumerate(mm_angles):
    #             an.value = scen[yidx]
    #         for s_idx, sc in enumerate(scan_scales):
    #             if args.bayes_eval_bin is not None:
    #                 if s_idx in args.bayes_eval_bin:
    #                     if idx == 0:
    #                         out += '_scale_{0:.0E}'.format(np.power(10, sc))
    #                 else: continue

    #             print '== SCALE = {0:.0E}'.format(np.power(10, sc))
    #             sc_angles.value = sc
    #             def lnProb(cube, ndim, nparams):
    #                 for i in range(ndim):
    #                     prange = prior_ranges[i][1] - prior_ranges[i][0]
    #                     p[i].value = prange*cube[i] + prior_ranges[i][0]
    #                 for name in p.names:
    #                     mcmc_paramset[name].value = p[name].value
    #                 theta = mcmc_paramset.values
    #                 # print 'theta', theta
    #                 # print 'mcmc_paramset', mcmc_paramset
    #                 return llh_utils.triangle_llh(
    #                     theta=theta,
    #                     args=args,
    #                     asimov_paramset=asimov_paramset,
    #                     mcmc_paramset=mcmc_paramset,
    #                     fitter=fitter
    #                 )
    #             prefix = 'mnrun/' + os.path.basename(out) + '_'
    #             if args.bayes_eval_bin is not None:
    #                 prefix += '{0}_{1}_'.format(idx, s_idx)
    #             print 'begin running evidence calculation for {0}'.format(prefix)
    #             result = pymultinest.run(
    #                 LogLikelihood=lnProb,
    #                 Prior=CubePrior,
    #                 n_dims=n_params,
    #                 importance_nested_sampling=True,
    #                 n_live_points=args.bayes_live_points,
    #                 evidence_tolerance=args.bayes_tolerance,
    #                 outputfiles_basename=prefix,
    #                 resume=False,
    #                 verbose=True
    #             )

    #             analyzer = pymultinest.Analyzer(outputfiles_basename=prefix, n_params=n_params)
    #             a_lnZ = analyzer.get_stats()['global evidence']
    #             print 'Evidence = {0}'.format(a_lnZ)
    #             scales.append(sc)
    #             evidences.append(a_lnZ)

    #         for i, d in enumerate(evidences):
    #             data[idx][i][0] = scales[i]
    #             data[idx][i][1] = d

    #     misc_utils.make_dir(out)
    #     print 'saving to {0}.npy'.format(out)
    #     np.save(out+'.npy', np.array(data))

    # dirname = os.path.dirname(out)
    # plot_utils.plot_BSM_angles_limit(
    #     dirname=dirname, outfile=outfile, outformat=['png'],
    #     args=args, bayesian=True
    # )

    # out = args.angles_corr_output+'/fr_co_evidence' + misc_utils.gen_identifier(args)
    # if args.run_angles_correlation:
    #     if args.bayes_eval_bin is None: assert 0
    #     import pymultinest

    #     scenarios = [
    #         [1, 0, 0, 0],
    #         [0, 1, 0, 0],
    #         [0, 0, 1, 0],
    #     ]
    #     nuisance = mcmc_paramset.from_tag(ParamTag.NUISANCE)
    #     mm_angles = mcmc_paramset.from_tag(ParamTag.MMANGLES)
    #     sc_angles = mcmc_paramset.from_tag(ParamTag.SCALE)[0]

    #     if not args.run_mcmc and args.likelihood is Likelihood.GOLEMFIT:
    #         fitter = gf_utils.setup_fitter(args, asimov_paramset)
    #     else: fitter = None

    #     def CubePrior(cube, ndim, nparams):
    #         # default are uniform priors
    #         return ;

    #     scan_angles = np.linspace(0, 1, args.bayes_bins)

    #     if args.bayes_eval_bin is not None:
    #         data = np.zeros((len(scenarios), 1, 1, 3))
    #     else: data = np.zeros((len(scenarios), args.bayes_bins, args.bayes_bins, 3))
    #     for idx, scen in enumerate(scenarios):
    #         for an in mm_angles:
    #             an.value = 0
    #         keep = mcmc_paramset.from_tag(ParamTag.MMANGLES)[idx]
    #         q = ParamSet(nuisance.params + [x for x in mm_angles if x.name != keep.name])
    #         n_params = len(q)
    #         prior_ranges = q.seeds

    #         scales, angles, evidences = [], [], []
    #         for s_idx, sc in enumerate(scan_scales):
    #             for a_idx, an in enumerate(scan_angles):
    #                 index = s_idx*args.bayes_bins + a_idx
    #                 if args.bayes_eval_bin is not None:
    #                     if index in args.bayes_eval_bin:
    #                         if idx == 0:
    #                             out += '_idx_{0}'.format(index)
    #                     else: continue

    #                 print '== SCALE = {0:.0E}'.format(np.power(10, sc))
    #                 print '== ANGLE = {0:.2f}'.format(an)
    #                 sc_angles.value = sc
    #                 keep.value = an
    #                 def lnProb(cube, ndim, nparams):
    #                     for i in range(ndim-1):
    #                         prange = prior_ranges[i][1] - prior_ranges[i][0]
    #                         q[i].value = prange*cube[i] + prior_ranges[i][0]
    #                     for name in q.names:
    #                         mcmc_paramset[name].value = q[name].value
    #                     theta = mcmc_paramset.values
    #                     # print 'theta', theta
    #                     # print 'mcmc_paramset', mcmc_paramset
    #                     return llh_utils.triangle_llh(
    #                         theta=theta,
    #                         args=args,
    #                         asimov_paramset=asimov_paramset,
    #                         mcmc_paramset=mcmc_paramset,
    #                         fitter=fitter
    #                     )
    #                 prefix = 'mnrun/' + os.path.basename(out) + '_'
    #                 if args.bayes_eval_bin is not None:
    #                     prefix += '{0}_{1}_{2}'.format(idx, s_idx, a_idx)

    #                 print 'begin running evidence calculation for {0}'.format(prefix)
    #                 result = pymultinest.run(
    #                     LogLikelihood=lnProb,
    #                     Prior=CubePrior,
    #                     n_dims=n_params,
    #                     importance_nested_sampling=True,
    #                     n_live_points=args.bayes_live_points,
    #                     evidence_tolerance=args.bayes_tolerance,
    #                     outputfiles_basename=prefix,
    #                     resume=False,
    #                     verbose=True
    #                 )

    #                 analyzer = pymultinest.Analyzer(outputfiles_basename=prefix, n_params=n_params)
    #                 a_lnZ = analyzer.get_stats()['global evidence']
    #                 print 'Evidence = {0}'.format(a_lnZ)
    #                 scales.append(sc)
    #                 angles.append(an)
    #                 evidences.append(a_lnZ)

    #             for i, d in enumerate(evidences):
    #                 data[idx][i][i][0] = scales[i]
    #                 data[idx][i][i][1] = angles[i]
    #                 data[idx][i][i][2] = d

    #     misc_utils.make_dir(out)
    #     print 'saving to {0}.npy'.format(out)
    #     np.save(out+'.npy', np.array(data))


main.__doc__ = __doc__


if __name__ == '__main__':
    main()

