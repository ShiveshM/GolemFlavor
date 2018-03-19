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
from collections import Sequence
import multiprocessing

import numpy as np

from utils.enums import Likelihood


class Param(object):
    """Parameter class to store parameters.
    """
    def __init__(self, name, value, ranges, std=None, tex=None):
        self._ranges = None
        self._tex = None

        self.name = name
        self.value = value
        self.ranges = ranges
        self.std = std
        self.tex = tex

    @property
    def ranges(self):
        return tuple(self._ranges)

    @ranges.setter
    def ranges(self, values):
        self._ranges = [val for val in values]

    @property
    def tex(self):
        return r'{0}'.format(self.tex)

    @tex.setter
    def tex(self, t):
        self._tex = t if t is not None else r'{\rm %s}' % self.name


class ParamSet(Sequence):
    """Container class for a set of parameters.
    """
    def __init__(self, *args):
        param_sequence = []
        for arg in args:
            try:
                param_sequence.extend(arg)
            except TypeError:
                param_sequence.append(arg)

        # Disallow duplicated params
        all_names = [p.name for p in param_sequence]
        unique_names = set(all_names)
        if len(unique_names) != len(all_names):
            duplicates = set([x for x in all_names if all_names.count(x) > 1])
            raise ValueError('Duplicate definitions found for param(s): ' +
                             ', '.join(str(e) for e in duplicates))

        # Elements of list must be Param type
        assert all([isinstance(x, Param) for x in param_sequence]), \
                'All params must be of type "Param"'

        self._params = param_sequence

    def __len__(self):
        return len(self._params)

    def __getitem__(self, i):
        if isinstance(i, int):
            return self._params[i]
        elif isinstance(i, basestring):
            return self._by_name[i]

    def __getattr__(self, attr):
        try:
            return super(ParamSet, self).__getattribute__(attr)
        except AttributeError:
            t, v, tb = sys.exc_info()
            try:
                return self[attr]
            except KeyError:
                raise t, v, tb

    def __iter__(self):
        return iter(self._params)

    @property
    def _by_name(self):
        return {obj.name: obj for obj in self._params}

    @property
    def names(self):
        return tuple([obj.name for obj in self._params])

    @property
    def values(self):
        return tuple([obj.value for obj in self._params])

    @property
    def ranges(self):
        return tuple([obj.ranges for obj in self._params])

    @property
    def stds(self):
        return tuple([obj.std for obj in self._params])

    @property
    def params(self):
        return self._params

    def to_dict(self):
        return {obj.name: obj.value for obj in self._params}


def gen_outfile_name(args):
    """Generate a name for the output file based on the input args.

    Parameters
    ----------
    args : argparse
        argparse object to print

    """
    mr = args.measured_ratio
    si = args.sigma_ratio
    if args.fix_source_ratio:
        sr = args.source_ratio
        if args.fix_mixing:
            outfile = args.outfile+'_{0:03d}_{1:03d}_{2:03d}_{3:04d}_sfr_{4:03d}_{5:03d}_{6:03d}_DIM{7}_fix_mixing'.format(
                int(mr[0]*100), int(mr[1]*100), int(mr[2]*100), int(si*1000),
                int(sr[0]*100), int(sr[1]*100), int(sr[2]*100), args.dimension
            )
        elif args.fix_scale:
            outfile = args.outfile+'_{0:03d}_{1:03d}_{2:03d}_{3:04d}_sfr_{4:03d}_{5:03d}_{6:03d}_DIM{7}_fixed_scale_{8}'.format(
                int(mr[0]*100), int(mr[1]*100), int(mr[2]*100), int(si*1000),
                int(sr[0]*100), int(sr[1]*100), int(sr[2]*100), args.dimension,
                args.scale
            )
        else:
            outfile = args.outfile+'_{0:03d}_{1:03d}_{2:03d}_{3:04d}_sfr_{4:03d}_{5:03d}_{6:03d}_DIM{7}_single_scale'.format(
                int(mr[0]*100), int(mr[1]*100), int(mr[2]*100), int(si*1000),
                int(sr[0]*100), int(sr[1]*100), int(sr[2]*100), args.dimension
            )
    else:
        if args.fix_mixing:
            outfile = args.outfile+'_{0:03d}_{1:03d}_{2:03d}_{3:04d}_DIM{4}_fix_mixing'.format(
                int(mr[0]*100), int(mr[1]*100), int(mr[2]*100),
                int(si*1000), args.dimension
            )
        elif args.fix_scale:
            outfile = args.outfile+'_{0:03d}_{1:03d}_{2:03d}_{3:04d}_DIM{4}_fixed_scale_{5}'.format(
                int(mr[0]*100), int(mr[1]*100), int(mr[2]*100),
                int(si*1000), args.dimension, args.scale
            )
        else:
            outfile = args.outfile+'_{0:03d}_{1:03d}_{2:03d}_{3:04d}_DIM{4}'.format(
                int(mr[0]*100), int(mr[1]*100), int(mr[2]*100),
                int(si*1000), args.dimension
            )
    if args.likelihood is Likelihood.FLAT: outfile += '_flat'
    return outfile


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
    for key in arg_vars.iterkeys():
        print '== {0:<25} = {1}'.format(key, arg_vars[key])


def enum_keys(e):
    return e.__members__.keys()


def enum_parse(s, c):
    return c[s.upper()]


def thread_type(t):
    if t.lower() == 'max':
        return multiprocessing.cpu_count()
    else:
        return int(t)
