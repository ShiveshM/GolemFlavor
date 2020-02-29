# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : March 17, 2018

"""
Misc functions for the BSM flavour ratio analysis
"""

from __future__ import absolute_import, division, print_function

import os
import errno
import multiprocessing
from fractions import gcd

import argparse
from operator import attrgetter

import numpy as np

from golemflavor.enums import str_enum
from golemflavor.enums import DataType, Likelihood, Texture


class SortingHelpFormatter(argparse.HelpFormatter):
    """Sort argparse help options alphabetically."""
    def add_arguments(self, actions):
        actions = sorted(actions, key=attrgetter('option_strings'))
        super(SortingHelpFormatter, self).add_arguments(actions)


def solve_ratio(fr):
    denominator = reduce(gcd, fr)
    f = [int(x/denominator) for x in fr]
    allow = (1, 2, 0)
    if f[0] not in allow or f[1] not in allow or f[2] not in allow:
        return '{0:.2f}_{1:.2f}_{2:.2f}'.format(fr[0], fr[1], fr[2])
    else:
        return '{0}_{1}_{2}'.format(f[0], f[1], f[2])


def gen_identifier(args):
    f = '_DIM{0}'.format(args.dimension)
    f += '_sfr_' + solve_ratio(args.source_ratio)
    if args.data in [DataType.ASIMOV, DataType.REALISATION]:
        f += '_mfr_' + solve_ratio(args.injected_ratio)
    if args.texture is not Texture.NONE:
        f += '_{0}'.format(str_enum(args.texture))
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
    >>> print(parse_bool('true'))
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
        print('== {0:<25} = {1}'.format(key, arg_vars[key]))


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


def centers(x):
    return (x[:-1]+x[1:])*0.5


def get_units(dimension):
    if dimension == 3: return r' / \:{\rm GeV}'
    if dimension == 4: return r''
    if dimension == 5: return r' / \:{\rm GeV}^{-1}'
    if dimension == 6: return r' / \:{\rm GeV}^{-2}'
    if dimension == 7: return r' / \:{\rm GeV}^{-3}'
    if dimension == 8: return r' / \:{\rm GeV}^{-4}'


def calc_nbins(x):
    n =  (np.max(x) - np.min(x)) / (2 * len(x)**(-1./3) * (np.percentile(x, 75) - np.percentile(x, 25)))
    return np.floor(n)


def calc_bins(x):
    nbins = calc_nbins(x)
    return np.linspace(np.min(x), np.max(x)+2, num=nbins+1)


def most_likely(arr):
    """Return the densest region given a 1D array of data."""
    binning = calc_bins(arr)
    harr = np.histogram(arr, binning)[0]
    return centers(binning)[np.argmax(harr)]


def interval(arr, percentile=68.):
    """Returns the *percentile* shortest interval around the mode."""
    center = most_likely(arr)
    sarr = sorted(arr)
    delta = np.abs(sarr - center)
    curr_low = np.argmin(delta)
    curr_up = curr_low
    npoints = len(sarr)
    while curr_up - curr_low < percentile/100.*npoints:
        if curr_low == 0:
            curr_up += 1
        elif curr_up == npoints-1:
            curr_low -= 1
        elif sarr[curr_up]-sarr[curr_low-1] < sarr[curr_up+1]-sarr[curr_low]:
            curr_low -= 1
        elif sarr[curr_up]-sarr[curr_low-1] > sarr[curr_up+1]-sarr[curr_low]:
            curr_up += 1
        elif (curr_up - curr_low) % 2:
            # they are equal so step half of the time up and down
            curr_low -= 1
        else:
            curr_up += 1
    return sarr[curr_low], center, sarr[curr_up]


def myround(x, base=2, up=False, down=False):
    if up == down and up is True: assert 0
    if up: return int(base * np.round(float(x)/base-0.5))
    elif down: return int(base * np.round(float(x)/base+0.5))
    else: int(base * np.round(float(x)/base))


