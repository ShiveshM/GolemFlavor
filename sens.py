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

import glob

import numpy as np
import numpy.ma as ma
from scipy.optimize import minimize

from utils import fr as fr_utils
from utils import gf as gf_utils
from utils import llh as llh_utils
from utils import misc as misc_utils
from utils import mn as mn_utils
from utils.enums import str_enum
from utils.enums import DataType, Likelihood, ParamTag
from utils.enums import PriorsCateg, StatCateg, Texture
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
    tag = ParamTag.NUISANCE
    nuisance.extend([
        Param(name='convNorm',        value=1.,  seed=[0.5, 2. ], ranges=[0.1, 10.], std=0.4, prior=lg_prior, tag=tag),
        Param(name='promptNorm',      value=0.,  seed=[0. , 6. ], ranges=[0. , 20.], std=2.4, prior=lg_prior, tag=tag),
        Param(name='muonNorm',        value=1.,  seed=[0.1, 2. ], ranges=[0. , 10.], std=0.1, tag=tag),
        Param(name='astroNorm',       value=8.0, seed=[0.,  5. ], ranges=[0. , 20.], std=1.5, tag=tag),
        Param(name='astroDeltaGamma', value=2.5, seed=[2.4, 3. ], ranges=[-5., 5. ], std=0.1, tag=tag)
    ])
    return ParamSet(nuisance)


def get_paramsets(args, nuisance_paramset):
    """Make the paramsets for generating the Asmimov MC sample and also running
    the MCMC.
    """
    asimov_paramset = []
    llh_paramset = []

    gf_nuisance = [x for x in nuisance_paramset.from_tag(ParamTag.NUISANCE)]

    llh_paramset.extend(
        [x for x in nuisance_paramset.from_tag(ParamTag.SM_ANGLES)]
    )
    llh_paramset.extend(gf_nuisance)

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
    if args.data is not DataType.REAL:
        flavour_angles = fr_utils.fr_to_angles(args.injected_ratio)
    else:
        flavour_angles = fr_utils.fr_to_angles([1, 1, 1])

    asimov_paramset.extend(gf_nuisance)
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
    args.source_ratio = fr_utils.normalise_fr(args.source_ratio)
    if args.data is not DataType.REAL:
        args.injected_ratio = fr_utils.normalise_fr(args.injected_ratio)

    args.binning = np.logspace(
        np.log10(args.binning[0]), np.log10(args.binning[1]), args.binning[2]+1
    )

    if args.eval_segment.lower() == 'all':
        args.eval_segment = None
    else:
        args.eval_segment = int(args.eval_segment)

    if args.stat_method is StatCateg.BAYESIAN:
        args.likelihood = Likelihood.GOLEMFIT
    elif args.stat_method is StatCateg.FREQUENTIST:
        args.likelihood = Likelihood.GF_FREQ

    if args.texture is Texture.NONE:
        raise ValueError('Must assume a BSM texture')


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
        '--datadir', type=str, default='./untitled',
        help='Path to store chains'
    )
    parser.add_argument(
        '--segments', type=int, default=10,
        help='Number of new physics scales to evaluate'
    )
    parser.add_argument(
        '--eval-segment', type=str, default='all',
        help='Which point to evalaute'
    )
    parser.add_argument(
        '--overwrite', type=misc_utils.parse_bool, default='False',
        help='Overwrite chains'
    )
    fr_utils.fr_argparse(parser)
    gf_utils.gf_argparse(parser)
    llh_utils.llh_argparse(parser)
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

    # Scale and BSM mixings will be fixed.
    scale_prm    = llh_paramset.from_tag(ParamTag.SCALE)[0]
    base_mn_pset = llh_paramset.from_tag(ParamTag.SCALE, invert=True)

    # Array of scales to scan over.
    boundaries = fr_utils.SCALE_BOUNDARIES[args.dimension]
    eval_scales = np.linspace(boundaries[0], boundaries[1], args.segments-1)
    eval_scales = np.concatenate([[-100.], eval_scales])

    # Evaluate just one point (job), or all points.
    if args.eval_segment is None:
        eval_dim = args.segments
    else: eval_dim = 1

    outfile = args.datadir + '/{0}/{1}/fr_stat'.format(
        *map(misc_utils.parse_enum, [args.stat_method, args.data])
    ) + misc_utils.gen_identifier(args)

    if not args.overwrite and os.path.isfile(outfile+'.npy'):
        print 'FILE EXISTS {0}'.format(outfile+'.npy')
        print 'Exiting...'
        return

    # Setup Golemfit.
    if args.run_mn:
        gf_utils.setup_fitter(args, asimov_paramset)

    # Initialise data structure.
    stat_arr = np.full((eval_dim, 2), np.nan)

    for idx_sc, scale in enumerate(eval_scales):
        if args.eval_segment is not None:
            if idx_sc == args.eval_segment:
                outfile += '_scale_{0:.0E}'.format(np.power(10, scale))
            else: continue
        print '|||| SCALE = {0:.0E}'.format(np.power(10, scale))

        if not args.overwrite and os.path.isfile(outfile+'.npy'):
            print 'FILE EXISTS {0}'.format(outfile+'.npy')
            t = np.load(outfile+'.npy')
            if np.any(~np.isfinite(t)):
                print 'nan found, rerunning...'
                pass
            else:
                print 'Exiting...'
                return

        # Lower scale boundary for first (NULL) point and set the scale param.
        reset_range = None
        if scale < scale_prm.ranges[0]:
            reset_range = scale_prm.ranges
            scale_prm.ranges = (scale, scale_prm.ranges[1])
        scale_prm.value = scale

        identifier = 'b{0}_{1}_{2}_sca{3}'.format(
            args.eval_segment, args.segments, str_enum(args.texture), scale
        )
        llh = '{0}'.format(args.likelihood).split('.')[1]
        data = '{0}'.format(args.data).split('.')[1]
        src_string = misc_utils.solve_ratio(args.source_ratio)
        prefix = args.mn_output + '/DIM{0}/{1}/{2}/s{3}/{4}'.format(
            args.dimension, data, llh, src_string, identifier
        )
        try:
            stat = mn_utils.mn_evidence(
                mn_paramset     = base_mn_pset,
                llh_paramset    = llh_paramset,
                asimov_paramset = asimov_paramset,
                args            = args,
                prefix          = prefix
            )
        except:
            print 'Failed run'
            raise
        print '## Evidence = {0}'.format(stat)

        if args.eval_segment is not None:
            stat_arr[0] = np.array([scale, stat])
        else:
            stat_arr[idx_sc] = np.array([scale, stat])

        # Cleanup.
        if reset_range is not None:
            scale_prm.ranges = reset_range

        if args.run_mn and not args.debug:
            try:
                for f in glob.glob(prefix + '*'):
                    print 'cleaning file {0}'.format(f)
                    os.remove(f)
            except:
                print 'got error trying to cleanup, continuing'
                pass

    misc_utils.make_dir(outfile)
    print 'Saving to {0}'.format(outfile+'.npy')
    np.save(outfile+'.npy', stat_arr)


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
