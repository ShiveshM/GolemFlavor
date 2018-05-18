#! /usr/bin/env python
# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : April 28, 2018

"""
HESE BSM flavour ratio analysis plotting script
"""

from __future__ import absolute_import, division

import os
import argparse
from functools import partial
from copy import deepcopy

import numpy as np
import numpy.ma as ma

from utils import fr as fr_utils
from utils import gf as gf_utils
from utils import likelihood as llh_utils
from utils import misc as misc_utils
from utils import plot as plot_utils
from utils.enums import EnergyDependance, Likelihood, ParamTag
from utils.enums import PriorsCateg, SensitivityCateg, StatCateg
from utils.param import Param, ParamSet, get_paramsets

from utils import multinest as mn_utils

import matplotlib as mpl
print mpl.rcParams.keys()
mpl.rcParams['text.usetex'] = True
# mpl.rcParams['text.latex.unicode'] = True
# mpl.rcParams['text.latex.preamble'] = r'\usepackage{cmbright}'
mpl.rcParams['font.family'] = 'sans-serif'
# mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
# mpl.rcParams['mathtext.fontset'] = 'stixsans'
# mpl.rcParams['mathtext.rm'] = 'DejaVu Sans'
# mpl.rcParams['mathtext.it'] = 'DejaVu Sans:italic'
# mpl.rcParams['mathtext.bf'] = 'DejaVu Sans:bold'


def define_nuisance():
    """Define the nuisance parameters."""
    tag = ParamTag.SM_ANGLES
    g_prior = PriorsCateg.GAUSSIAN
    hg_prior = PriorsCateg.HALFGAUSS
    e = 1e-9
    nuisance = [
        Param(name='s_12_2', value=0.307,          seed=[0.26, 0.35],     ranges=[0., 1.],      std=0.013,   tex=r's_{12}^2', prior=g_prior,  tag=tag),
        Param(name='c_13_4', value=1-(0.02206)**2, seed=[0.995, 1-e],     ranges=[0., 1.],      std=0.00147, tex=r'c_{13}^4', prior=hg_prior, tag=tag),
        Param(name='s_23_2', value=0.538,          seed=[0.31, 0.75],     ranges=[0., 1.],      std=0.069,   tex=r's_{23}^2', prior=g_prior,  tag=tag),
        Param(name='dcp',    value=4.08404,        seed=[0+e, 2*np.pi-e], ranges=[0., 2*np.pi], std=2.0,     tex=r'\delta_{CP}', tag=tag),
        Param(
            name='m21_2', value=7.40E-23, seed=[7.2E-23, 7.6E-23], ranges=[6.80E-23, 8.02E-23],
            std=2.1E-24, tex=r'\Delta m_{21}^2{\rm GeV}^{-2}', prior=g_prior, tag=tag
        ),
        Param(
            name='m3x_2', value=2.494E-21, seed=[2.46E-21, 2.53E-21], ranges=[2.399E-21, 2.593E-21],
            std=3.3E-23, tex=r'\Delta m_{3x}^2{\rm GeV}^{-2}', prior=g_prior, tag=tag
        )
    ]
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
    if args.fix_mixing and args.fix_scale:
        raise NotImplementedError('Fixed mixing and scale not implemented')
    if args.fix_mixing and args.fix_mixing_almost:
        raise NotImplementedError(
            '--fix-mixing and --fix-mixing-almost cannot be used together'
        )
    if args.fix_scale:
        raise NotImplementedError(
            '--fix-scale not implemented'
        )

    args.measured_ratio = fr_utils.normalise_fr(args.measured_ratio)
    if args.fix_source_ratio:
        assert len(args.source_ratios) % 3 == 0
        srs = [args.source_ratios[3*x:3*x+3]
               for x in range(int(len(args.source_ratios)/3))]
        args.source_ratios = map(fr_utils.normalise_fr, srs)

    if args.energy_dependance is EnergyDependance.SPECTRAL:
        args.binning = np.logspace(
            np.log10(args.binning[0]), np.log10(args.binning[1]), args.binning[2]+1
        )

    if args.split_jobs and args.run_method is SensitivityCateg.FULL:
        raise NotImplementedError(
            'split_jobs and run_method not implemented'
        )

    args.dimensions = np.sort(args.dimensions)

    args_copy = deepcopy(args)
    scale_regions = []
    for dim in args.dimensions:
        args_copy.dimension = dim
        _, scale_region = fr_utils.estimate_scale(args_copy)
        scale_regions.append(scale_region)
    args.scale_region = [np.min(scale_regions), np.max(scale_regions)]
    args.scale = np.power(10., np.average(np.log10(args.scale_region)))


def parse_args(args=None):
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="HESE BSM flavour ratio analysis plotting script",
        formatter_class=misc_utils.SortingHelpFormatter,
    )
    parser.add_argument(
        '--infile', type=str, default='./untitled',
        help='Path to input dir'
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
        '--split-jobs', type=misc_utils.parse_bool, default='False',
        help='Did the jobs get split'
    )
    parser.add_argument(
        '--plot', type=misc_utils.parse_bool, default='True',
        help='Make sensitivity plots'
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
    misc_utils.remove_option(parser, 'dimension')
    misc_utils.remove_option(parser, 'source_ratio')
    misc_utils.remove_option(parser, 'scale')
    misc_utils.remove_option(parser, 'scale_region')
    parser.add_argument(
        '--dimensions', type=int, nargs='*', default=[3, 6],
        help='Set the new physics dimensions to consider'
    )
    parser.add_argument(
        '--source-ratios', type=int, nargs='*', default=[2, 1, 0],
        help='Set the source flavour ratios for the case when you want to fix it'
    )
    if args is None: return parser.parse_args()
    else: return parser.parse_args(args.split())


def main():
    args = parse_args()
    process_args(args)
    misc_utils.print_args(args)

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

    corr_angles_categ = [SensitivityCateg.CORR_ANGLE, SensitivityCateg.CORR_ONE_ANGLE]
    fixed_angle_categ = [SensitivityCateg.FIXED_ANGLE, SensitivityCateg.FIXED_ONE_ANGLE]

    if args.run_method in corr_angles_categ:
        scan_angles = np.linspace(0+1e-9, 1-1e-9, args.sens_bins)
    else: scan_angles = np.array([0])
    print 'scan_angles', scan_angles

    dims = len(args.dimensions)
    srcs = len(args.source_ratios)
    if args.run_method is SensitivityCateg.FULL:
        statistic_arr = np.full((dims, srcs, args.sens_bins, 2), np.nan)
    elif args.run_method in fixed_angle_categ:
        statistic_arr = np.full((dims, srcs, len(st_paramset_arr), args.sens_bins, 2), np.nan)
    elif args.run_method in corr_angles_categ:
        statistic_arr = np.full(
            (dims, srcs, len(st_paramset_arr), args.sens_bins, args.sens_bins, 3), np.nan
        )

    print 'Loading data'
    for idim, dim in enumerate(args.dimensions):
        argsc = deepcopy(args)
        argsc.dimension = dim
        _, scale_region = fr_utils.estimate_scale(argsc)
        argsc.scale_region = scale_region
        scan_scales = np.linspace(
            np.log10(scale_region[0]), np.log10(scale_region[1]), args.sens_bins
        )

        for isrc, src in enumerate(args.source_ratios):
            argsc.source_ratio = src
            infile = args.infile
            if args.likelihood is Likelihood.GOLEMFIT:
                infile += '/golemfit/'
            elif args.likelihood is Likelihood.GAUSSIAN:
                infile += '/gaussian/'
            if args.likelihood is Likelihood.GAUSSIAN:
                infile += '{0}/'.format(str(args.sigma_ratio).replace('.', '_'))
            infile += '/DIM{0}/fix_ifr/{1}/{2}/{3}/fr_stat'.format(
            # infile += '/DIM{0}/fix_ifr/100TeV/{1}/{2}/{3}/fr_stat'.format(
                dim, *map(misc_utils.parse_enum, [args.stat_method, args.run_method, args.data])
            ) + misc_utils.gen_identifier(argsc)
            print '== {0:<25} = {1}'.format('infile', infile)

            if args.split_jobs:
                for idx_an, an in enumerate(scan_angles):
                    for idx_sc, sc in enumerate(scan_scales):
                        filename = infile + '_scale_{0:.0E}'.format(np.power(10, sc))
                        try:
                            if args.run_method in fixed_angle_categ:
                                print 'Loading from {0}'.format(filename+'.npy')
                                statistic_arr[idim][isrc][:,idx_sc] = np.load(filename+'.npy')[:,0]
                            if args.run_method in corr_angles_categ:
                                filename += '_angle_{0:<04.2}'.format(an)
                                print 'Loading from {0}'.format(filename+'.npy')
                                statistic_arr[idim][isrc][:,idx_an,idx_sc] = np.load(filename+'.npy')[:,0,0]
                        except:
                            print 'Unable to load file {0}'.format(filename+'.npy')
                            continue
            else:
                print 'Loading from {0}'.format(infile+'.npy')
                try:
                    statistic_arr[idim][isrc] = np.load(infile+'.npy')
                except:
                    print 'Unable to load file {0}'.format(infile+'.npy')
                    continue

    print 'statistic_arr', statistic_arr

    data = ma.masked_invalid(statistic_arr)
    if args.plot_statistic:
        print 'Plotting statistic'

        argsc = deepcopy(args)
        for idim, dim in enumerate(args.dimensions):
            argsc.dimension = dim
            _, scale_region = fr_utils.estimate_scale(argsc)
            argsc.scale_region = scale_region
            base_infile = args.infile
            if args.likelihood is Likelihood.GOLEMFIT:
                base_infile += '/golemfit/'
            elif args.likelihood is Likelihood.GAUSSIAN:
                base_infile += '/gaussian/'
            if args.likelihood is Likelihood.GAUSSIAN:
                base_infile += '{0}/'.format(str(args.sigma_ratio).replace('.', '_'))
            base_infile += '/DIM{0}/fix_ifr'.format(dim)
            # base_infile += '/DIM{0}/fix_ifr/100TeV'.format(dim)

            for isrc, src in enumerate(args.source_ratios):
                argsc.source_ratio = src
                infile = base_infile +'/{0}/{1}/{2}/fr_stat'.format(
                    *map(misc_utils.parse_enum, [args.stat_method, args.run_method, args.data])
                ) + misc_utils.gen_identifier(argsc)
                basename = os.path.dirname(infile)
                baseoutfile = basename[:5]+basename[5:].replace('data', 'plots')
                baseoutfile += '/' + os.path.basename(infile)
                if args.run_method is SensitivityCateg.FULL:
                    outfile = baseoutfile
                    plot_utils.plot_statistic(
                        data        = data[idim][isrc],
                        outfile     = outfile,
                        outformat   = ['png'],
                        args        = argsc,
                        scale_param = scale,
                    )
                if args.run_method in fixed_angle_categ:
                    for idx_scen in xrange(len(st_paramset_arr)):
                        print '|||| SCENARIO = {0}'.format(idx_scen)
                        outfile = baseoutfile + '_SCEN{0}'.format(idx_scen)
                        if idx_scen == 0: label = r'$\mathcal{O}_{12}=\frac{\pi}{4}$'
                        elif idx_scen == 1: label = r'$\mathcal{O}_{13}=\frac{\pi}{4}$'
                        elif idx_scen == 2: label = r'$\mathcal{O}_{23}=\frac{\pi}{4}$'
                        plot_utils.plot_statistic(
                            data        = data[idim][isrc][idx_scen],
                            outfile     = outfile,
                            outformat   = ['png'],
                            args        = argsc,
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
                                data        = data[idim][isrc][idx_scen][idx_an][:,1:],
                                outfile     = outfile,
                                outformat   = ['png'],
                                args        = argsc,
                                scale_param = scale,
                                label       = _label
                            )

    if args.plot:
        print 'Plotting sensitivities'

        basename = args.infile[:5]+args.infile[5:].replace('data', 'plots')
        baseoutfile = basename + '/{0}/{1}/{2}/'.format(
            *map(misc_utils.parse_enum, [args.likelihood, args.stat_method, args.data])
        )

        if args.run_method is SensitivityCateg.FULL:
            plot_utils.plot_sens_full(
                data      = data,
                outfile   = baseoutfile + '/FULL',
                outformat = ['png'],
                args      = args,
            )
        elif args.run_method in fixed_angle_categ:
            plot_utils.plot_sens_fixed_angle(
                data      = data,
                outfile   = baseoutfile + '/FIXED_ANGLE',
                outformat = ['png'],
                args      = args,
            )
        elif args.run_method in corr_angles_categ:
            plot_utils.plot_sens_corr_angle(
                data      = data,
                outfile   = baseoutfile + '/CORR_ANGLE',
                outformat = ['png'],
                args      = args,
            )


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
