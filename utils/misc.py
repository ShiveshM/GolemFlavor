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
from fractions import gcd

import argparse
from operator import attrgetter

import numpy as np

from utils.enums import Likelihood


class SortingHelpFormatter(argparse.HelpFormatter):
    """Sort argparse help options alphabetically."""
    def add_arguments(self, actions):
        actions = sorted(actions, key=attrgetter('option_strings'))
        super(SortingHelpFormatter, self).add_arguments(actions)


def solve_ratio(fr):
    denominator = reduce(gcd, fr)
    return [int(x/denominator) for x in fr]


def gen_identifier(args):
    f = '_DIM{0}'.format(args.dimension)
    mr1, mr2, mr3 = solve_ratio(args.measured_ratio)
    if args.fix_source_ratio:
        sr1, sr2, sr3 = solve_ratio(args.source_ratio)
        f += '_sfr_{0:G}_{1:G}_{2:G}_mfr_{3:G}_{4:G}_{5:G}'.format(
            sr1, sr2, sr3, mr1, mr2, mr3
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


def parse_enum(e):
    return '{0}'.format(e).split('.')[1].lower()


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

def remove_option(parser, arg):
    for action in parser._actions:
        if (vars(action)['option_strings']
            and vars(action)['option_strings'][0] == arg) \
                or vars(action)['dest'] == arg:
            parser._remove_action(action)

    for action in parser._action_groups:
        vars_action = vars(action)
        var_group_actions = vars_action['_group_actions']
        for x in var_group_actions:
            if x.dest == arg:
                var_group_actions.remove(x)
                return


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

