#! /usr/bin/env python
"""
Sample points for a specific scenario
"""

from __future__ import absolute_import, division

import sys
sys.path.extend(['.', '../'])

import argparse
from functools import partial

import numpy as np

from utils import fr as fr_utils
from utils import misc as misc_utils
from utils.param import Param, ParamSet, get_paramsets


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
        '--plot-statistic', type=misc_utils.parse_bool, default='False',
        help='Plot MultiNest evidence or LLH value'
    )
    fr_utils.fr_argparse(parser)
    if args is None: return parser.parse_args()
    else: return parser.parse_args(args.split())


def main():
    args = parse_args()
    process_args(args)
    misc_utils.print_args(args)

    if args.seed is not None:
        np.random.seed(args.seed)

    asimov_paramset, llh_paramset = get_paramsets(args, ParamSet())
    outfile = misc_utils.gen_outfile_name(args)
    print '== {0:<25} = {1}'.format('outfile', outfile)


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
