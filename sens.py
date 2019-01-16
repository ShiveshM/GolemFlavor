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
from scipy.optimize import minimize

from utils import fr as fr_utils
from utils import gf as gf_utils
from utils import llh as llh_utils
from utils import misc as misc_utils
from utils import plot as plot_utils
from utils.enums import EnergyDependance, Likelihood, MixingScenario, ParamTag
from utils.enums import PriorsCateg, SensitivityCateg, StatCateg
from utils.param import Param, ParamSet, get_paramsets

from utils import mn as mn_utils


def define_nuisance():
    """Define the nuisance parameters."""
    tag = ParamTag.SM_ANGLES
    nuisance = []
    g_prior = PriorsCateg.GAUSSIAN
    e = 1e-9
    nuisance.extend([
        Param(name='s_12_2', value=0.307,            seed=[0.26, 0.35],     ranges=[0., 1.],      std=0.013,   tex=r's_{12}^2', prior=g_prior,  tag=tag),
        Param(name='c_13_4', value=(1-(0.02206))**2, seed=[0.950, 0.961],   ranges=[0., 1.],      std=0.00147, tex=r'c_{13}^4', prior=g_prior,  tag=tag),
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
        Param(name='astroDeltaGamma', value=2.5, seed=[2.4, 3. ], ranges=[-5., 5. ], std=0.1,  tag=tag)
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
    if args.sens_run and args.fix_scale:
        raise NotImplementedError(
            '--sens-run and --fix-scale cannot be used together'
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

    if args.sens_eval_bin.lower() == 'all':
        args.sens_eval_bin = None
    else:
        args.sens_eval_bin = int(args.sens_eval_bin)

    if args.sens_eval_bin is not None and args.plot_statistic:
        print 'Cannot make plot when specific scale bin is chosen'
        args.plot_statistic = False

    if args.stat_method is StatCateg.FREQUENTIST and \
       args.likelihood is Likelihood.GOLEMFIT:
        args.likelihood = Likelihood.GF_FREQ

    args.burnin = False


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
        '--sens-run', type=misc_utils.parse_bool, default='True',
        help='Generate sensitivities'
    )
    parser.add_argument(
        '--run-method', default='full',
        type=partial(misc_utils.enum_parse, c=SensitivityCateg),
        choices=SensitivityCateg,
        help='Choose which type of sensivity plot to make'
    )
    parser.add_argument(
        '--stat-method', default='bayesian',
        type=partial(misc_utils.enum_parse, c=StatCateg), choices=StatCateg,
        help='Statistical method to employ'
    )
    parser.add_argument(
        '--sens-bins', type=int, default=10,
        help='Number of bins for the Bayes factor plot'
    )
    parser.add_argument(
        '--sens-eval-bin', type=str, default='all',
        help='Which bin to evalaute for Bayes factor plot'
    )
    parser.add_argument(
        '--plot-statistic', type=misc_utils.parse_bool, default='False',
        help='Plot MultiNest evidence or LLH value'
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

    asimov_paramset, llh_paramset = get_paramsets(args, define_nuisance())

    scale = llh_paramset.from_tag(ParamTag.SCALE)[0]
    mmangles = llh_paramset.from_tag(ParamTag.MMANGLES)
    if args.run_method is SensitivityCateg.FULL:
        st_paramset_arr = [llh_paramset.from_tag(ParamTag.SCALE, invert=True)]
    elif args.run_method in [SensitivityCateg.FIXED_ANGLE, SensitivityCateg.CORR_ANGLE]:
        nscale_pmset = llh_paramset.from_tag([ParamTag.SCALE, ParamTag.MMANGLES], invert=True)
        st_paramset_arr = [nscale_pmset] * 3
    elif args.run_method in [SensitivityCateg.FIXED_ONE_ANGLE, SensitivityCateg.CORR_ONE_ANGLE]:
        nscale_pmset = llh_paramset.from_tag(ParamTag.SCALE, invert=True)
        st_paramset_arr = []
        for x in xrange(3):
            st_paramset_arr.append(
                ParamSet([prms for prms in nscale_pmset
                          if mmangles[x].name != prms.name])
            )

    scan_scales = np.linspace(
        np.log10(args.scale_region[0]), np.log10(args.scale_region[1]), args.sens_bins
    )
    scan_scales = np.concatenate([[-100.], scan_scales])
    # scan_scales = np.array([-100.])

    corr_angles_categ = [SensitivityCateg.CORR_ANGLE, SensitivityCateg.CORR_ONE_ANGLE]
    fixed_angle_categ = [SensitivityCateg.FIXED_ANGLE, SensitivityCateg.FIXED_ONE_ANGLE]

    if args.run_method in corr_angles_categ:
        scan_angles = np.linspace(0+1e-9, 1-1e-9, args.sens_bins)
    else: scan_angles = np.array([0])
    print 'scan_scales', scan_scales
    print 'scan_angles', scan_angles

    if args.sens_eval_bin is None:
        eval_dim = args.sens_bins
    else: eval_dim = 1

    out = args.outfile+'/{0}/{1}/{2}/fr_stat'.format(
        *map(misc_utils.parse_enum, [args.stat_method, args.run_method, args.data])
    ) + misc_utils.gen_identifier(args)
    if args.sens_run:
        if args.likelihood in [Likelihood.GOLEMFIT, Likelihood.GF_FREQ]:
            fitter = gf_utils.setup_fitter(args, asimov_paramset)
            if args.stat_method is StatCateg.FREQUENTIST:
                flags, gf_nuisance = gf_utils.fit_flags(llh_paramset)
                llh_paramset = llh_paramset.remove_params(gf_nuisance)
                asimov_paramset = asimov_paramset.remove_params(gf_nuisance)
                st_paramset_arr = [x.remove_params(gf_nuisance)
                                   for x in st_paramset_arr]
                fitter.SetFitParametersFlag(flags)
        else: fitter = None

        if args.run_method is SensitivityCateg.FULL:
            statistic_arr = np.full((eval_dim, 2), np.nan)
        elif args.run_method in fixed_angle_categ:
            statistic_arr = np.full((len(st_paramset_arr), eval_dim, 2), np.nan)
        elif args.run_method in corr_angles_categ:
            statistic_arr = np.full(
                (len(st_paramset_arr), eval_dim, eval_dim, 3), np.nan
            )

        for idx_scen, sens_paramset in enumerate(st_paramset_arr):
            print '|||| SCENARIO = {0}'.format(idx_scen)
            if args.run_method in fixed_angle_categ:
                for x in mmangles: x.value = 0.+1e-9
                if idx_scen == 0 or idx_scen == 2:
                    mmangles[idx_scen].value = np.sin(np.pi/4.)**2
                    """s_12^2 or s_23^2"""
                    mmangles[1].value = 1.
                    """c_13^4"""
                elif idx_scen == 1:
                    mmangles[idx_scen].value = np.cos(np.pi/4.)**4
                    """c_13^4"""

            for idx_an, an in enumerate(scan_angles):
                if args.run_method in corr_angles_categ:
                    for x in mmangles: x.value = 0.+1e-9
                    angle = np.arcsin(np.sqrt(an))/2.
                    if idx_scen == 0 or idx_scen == 2:
                        mmangles[idx_scen].value = np.sin(angle)**2
                        """s_12^2 or s_23^2"""
                        mmangles[1].value = 1.
                        """c_13^4"""
                    elif idx_scen == 1:
                        mmangles[idx_scen].value = np.cos(angle)**4
                        """c_13^4"""

                for idx_sc, sc in enumerate(scan_scales):
                    if args.sens_eval_bin is not None:
                        if args.run_method in corr_angles_categ:
                            index = idx_an*args.sens_bins + idx_sc
                        else: index = idx_sc
                        if index == args.sens_eval_bin:
                            if idx_scen == 0:
                                out += '_scale_{0:.0E}'.format(np.power(10, sc))
                                if args.run_method in corr_angles_categ:
                                    out += '_angle_{0:<04.2}'.format(an)
                        else: continue

                    sc_param = llh_paramset.from_tag(ParamTag.SCALE)[0]
                    if sc < sc_param.ranges[0]:
                        sc_param.ranges = (sc, sc_param.ranges[1])

                    if idx_sc == 0 or args.sens_eval_bin is not None:
                        print '|||| ANGLE = {0:<04.2}'.format(float(an))
                    print '|||| SCALE = {0:.0E}'.format(np.power(10, sc))
                    scale.value = sc
                    if args.stat_method is StatCateg.BAYESIAN:
                        identifier = 'b{0}_{1}_sce{2}_sca{3}_an{4}'.format(
                            args.sens_eval_bin, args.sens_bins, idx_scen, sc, idx_an
                        )
                        try:
                            stat = mn_utils.mn_evidence(
                                mn_paramset     = sens_paramset,
                                llh_paramset    = llh_paramset,
                                asimov_paramset = asimov_paramset,
                                args            = args,
                                fitter          = fitter,
                                identifier      = identifier
                            )
                        except:
                            print 'Failed run, continuing'
                            # raise
                            continue
                        print '## Evidence = {0}'.format(stat)
                    elif args.stat_method is StatCateg.FREQUENTIST:
                        def fn(x):
                            # TODO(shivesh): should be seed or ranges?
                            # Force prior ranges to be inside "seed"
                            for el in x:
                                if el < 0 or el > 1: return np.inf
                            pranges = sens_paramset.seeds
                            for i, name in enumerate(sens_paramset.names):
                                llh_paramset[name].value = \
                                    (pranges[i][1]-pranges[i][0])*x[i] + pranges[i][0]
                            theta = llh_paramset.values
                            llh = llh_utils.ln_prob(
                                theta=theta, args=args, asimov_paramset=asimov_paramset,
                                llh_paramset=llh_paramset, fitter=fitter
                            )
                            # print 'llh_paramset', llh_paramset
                            # print '-llh', -llh
                            return -llh

                        n_params = len(sens_paramset)
                        bounds = np.dstack([np.zeros(n_params), np.ones(n_params)])[0]
                        x0_arr = np.linspace(0+1e-3, 1-1e-3, 3)
                        stat = -np.inf
                        try:
                            for x0_v in x0_arr:
                                x0 = np.full(n_params, x0_v)
                                # res = minimize(fun=fn, x0=x0, method='Powell', tol=1e-3)
                                res = minimize(fun=fn, x0=x0, method='Nelder-Mead', tol=1e-3)
                                # res = minimize(fun=fn, x0=x0, method='L-BFGS-B', tol=1e-3, bounds=bounds)
                                # res = minimize(fun=fn, x0=x0, method='BFGS', tol=1e-3)
                                s = -fn(res.x)
                                if s > stat:
                                    stat = s
                        except AssertionError:
                            # print 'Failed run, continuing'
                            raise
                            continue
                        print '=== final llh', stat
                    if args.run_method is SensitivityCateg.FULL:
                        statistic_arr[idx_sc] = np.array([sc, stat])
                    elif args.run_method in fixed_angle_categ:
                        if args.sens_eval_bin is not None:
                            statistic_arr[idx_scen][0] = np.array([sc, stat])
                        else:
                            statistic_arr[idx_scen][idx_sc] = np.array([sc, stat])
                    elif args.run_method in corr_angles_categ:
                        if args.sens_eval_bin is not None:
                            statistic_arr[idx_scen][0][0] = np.array([an, sc, stat])
                        else:
                            statistic_arr[idx_scen][idx_an][idx_sc] = np.array([an, sc, stat])

        misc_utils.make_dir(out)
        print 'Saving to {0}'.format(out+'.npy')
        np.save(out+'.npy', statistic_arr)

    if args.plot_statistic:
        print 'Plotting statistic'
        if args.sens_run: raw = statistic_arr
        else: raw = np.load(out+'.npy')
        data = ma.masked_invalid(raw)

        basename = os.path.dirname(out) + '/statrun/' + os.path.basename(out)
        baseoutfile = basename[:5]+basename[5:].replace('data', 'plots')
        if args.run_method is SensitivityCateg.FULL:
            plot_utils.plot_statistic(
                data        = data,
                outfile     = baseoutfile,
                outformat   = ['png'],
                args        = args,
                scale_param = scale
            )
        elif args.run_method in fixed_angle_categ:
            for idx_scen in xrange(len(st_paramset_arr)):
                print '|||| SCENARIO = {0}'.format(idx_scen)
                outfile = baseoutfile + '_SCEN{0}'.format(idx_scen)
                if idx_scen == 0: label = r'$\mathcal{O}_{12}=\pi/4$'
                elif idx_scen == 1: label = r'$\mathcal{O}_{13}=\pi/4$'
                elif idx_scen == 2: label = r'$\mathcal{O}_{23}=\pi/4$'
                plot_utils.plot_statistic(
                    data        = data[idx_scen],
                    outfile     = outfile,
                    outformat   = ['png'],
                    args        = args,
                    scale_param = scale,
                    label       = label
                )
        elif args.run_method in corr_angles_categ:
            for idx_scen in xrange(len(st_paramset_arr)):
                print '|||| SCENARIO = {0}'.format(idx_scen)
                basescenoutfile = baseoutfile + '_SCEN{0}'.format(idx_scen)
                if idx_scen == 0: label = r'$\mathcal{O}_{12}='
                elif idx_scen == 1: label = r'$\mathcal{O}_{13}='
                elif idx_scen == 2: label = r'$\mathcal{O}_{23}='
                for idx_an, an in enumerate(scan_angles):
                    print '|||| ANGLE = {0:<04.2}'.format(float(an))
                    outfile = basescenoutfile + '_ANGLE{0}'.format(idx_an)
                    _label = label + r'{0:<04.2}$'.format(an)
                    plot_utils.plot_statistic(
                        data        = data[idx_scen][index][:,1:],
                        outfile     = outfile,
                        outformat   = ['png'],
                        args        = args,
                        scale_param = scale,
                        label       = _label
                    )


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
