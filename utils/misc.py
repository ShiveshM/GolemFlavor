# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : March 17, 2018

"""
Misc functions for the BSM flavour ratio analysis
"""

from __future__ import absolute_import, division

import os
import errno
import multiprocessing

import argparse
from operator import attrgetter

import numpy as np

from utils.enums import Likelihood


DTYPE  = np.float128
CDTYPE = np.complex128
PI     = np.arccos(DTYPE(-1))


class SortingHelpFormatter(argparse.HelpFormatter):
    """Sort argparse help options alphabetically."""
    def add_arguments(self, actions):
        actions = sorted(actions, key=attrgetter('option_strings'))
        super(SortingHelpFormatter, self).add_arguments(actions)


def gen_identifier(args):
    f = '_DIM{0}'.format(args.dimension)
    mr1, mr2, mr3 = args.measured_ratio
    if args.fix_source_ratio:
        sr1, sr2, sr3 = args.source_ratio
        f += '_sfr_{0:03d}_{1:03d}_{2:03d}_mfr_{3:03d}_{4:03d}_{5:03d}'.format(
            int(sr1*100), int(sr2*100), int(sr3*100),
            int(mr1*100), int(mr2*100), int(mr3*100)
        )
        if args.fix_mixing:
            f += '_fix_mixing'
        elif args.fix_mixing_almost:
            f += '_fix_mixing_almost'
        elif args.fix_scale:
            f += '_fix_scale_{0}'.format(args.scale)
    else:
        f += '_mfr_{3:03d}_{4:03d}_{5:03d}'.format(mr1, mr2, mr3)
    if args.fix_mixing:
        f += '_fix_mixing'
    elif args.fix_mixing_almost:
        f += '_fix_mixing_almost'
    elif args.fix_scale:
        f += '_fix_scale_{0}'.format(args.scale)
    if args.likelihood is Likelihood.FLAT: f += '_flat'
    elif args.likelihood is Likelihood.GAUSSIAN:
        f += '_sigma_{0:03d}'.format(int(args.sigma_ratio*1000))
    return f


def gen_outfile_name(args):
    """Generate a name for the output file based on the input args.

    Parameters
    ----------
    args : argparse
        argparse object to print

    """
    return args.outfile + gen_identifier(args)


def parse_bool(s):
    """Parse a string to a boolean.

    Parameters
    ----------
    s : str
        String to parse

    Returns
    ----------
    bool

    Examples
    ----------
    >>> from misc import parse_bool
    >>> print parse_bool('true')
    True

    """
    if s.lower() == 'true':
        return True
    elif s.lower() == 'false':
        return False
    else:
        raise ValueError


def print_args(args):
    """Print the input arguments.

    Parameters
    ----------
    args : argparse
        argparse object to print

    """
    arg_vars = vars(args)
    for key in sorted(arg_vars):
        print '== {0:<25} = {1}'.format(key, arg_vars[key])


def enum_parse(s, c):
    return c[s.upper()]


def make_dir(outfile):
    try:
        os.makedirs(outfile[:-len(os.path.basename(outfile))])
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(outfile[:-len(os.path.basename(outfile))]):
            pass
        else:
            raise


def seed_parse(s):
    if s.lower() == 'none':
        return None
    else:
        return int(s)


def thread_type(t):
    if t.lower() == 'max':
        return multiprocessing.cpu_count()
    else:
        return int(t)


def thread_factors(t):
    for x in reversed(range(int(np.ceil(np.sqrt(t)))+1)):
        if t%x == 0:
            return (x, int(t/x))

