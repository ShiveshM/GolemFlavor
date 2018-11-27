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

from utils import gf as gf_utils
from utils import llh as llh_utils
from utils import misc as misc_utils
from utils import plot as plot_utils
from utils.enums import str_enum
from utils.enums import DataType, Likelihood, ParamTag
from utils.param import Param, ParamSet, get_paramsets


def define_nuisance():
    """Define the nuisance parameters."""
    nuisance = []
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
        '--run-scan', type=misc_utils.parse_bool, default='True',
        help='Do the scan from scratch'
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
    nuisance_argparse(parser)
    misc_utils.remove_option(parser, 'sigma_ratio')
    if args is None: return parser.parse_args()
    else: return parser.parse_args(args.split())


def gen_identifier(args):
    f = '_{0}_{1}'.format(*map(str_enum, (args.likelihood, args.data)))
    if args.data is not DataType.REAL:
        ir1, ir2, ir3 = solve_ratio(args.injected_ratio)
        f += '_INJ_{1:03d}_{2:03d}_{3:03d}'.format(ir1, ir2, ir3)
    return f


def main():
    args = parse_args()
    process_args(args)
    misc_utils.print_args(args)

    if args.seed is not None:
        np.random.seed(args.seed)

    asimov_paramset, llh_paramset = get_paramsets(args, define_nuisance())
    outfile = args.outfile + gen_identifier(args)
    print '== {0:<25} = {1}'.format('outfile', outfile)

    if args.run_scan:
        fitter = gf_utils.setup_fitter(args, asimov_paramset)

    print "DONE!"


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
