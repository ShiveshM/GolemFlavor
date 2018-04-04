# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : March 17, 2018

"""
Misc functions for the BSM flavour ratio analysis
"""

from __future__ import absolute_import, division

import os, sys
import errno
import multiprocessing

import argparse
from collections import Sequence
from operator import attrgetter

import numpy as np

from utils.enums import Likelihood, ParamTag


class Param(object):
    """Parameter class to store parameters.
    """
    def __init__(self, name, value, ranges, std=None, tex=None, tag=None):
        self._ranges = None
        self._tex = None
        self._tag = None

        self.name = name
        self.value = value
        self.ranges = ranges
        self.std = std
        self.tex = tex
        self.tag = tag

    @property
    def ranges(self):
        return tuple(self._ranges)

    @ranges.setter
    def ranges(self, values):
        self._ranges = [val for val in values]

    @property
    def tex(self):
        return r'{0}'.format(self._tex)

    @tex.setter
    def tex(self, t):
        self._tex = t if t is not None else r'{\rm %s}' % self.name

    @property
    def tag(self):
        return self._tag

    @tag.setter
    def tag(self, t):
        if t is None: self._tag = ParamTag.NONE
        else:
            assert t in ParamTag
            self._tag = t


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

    def __str__(self):
        o = '\n'
        for obj in self._params:
            o += '== {0:<15} = {1:<15}, tag={2:<15}\n'.format(
                obj.name, obj.value, obj.tag
            )
        return o

    @property
    def _by_name(self):
        return {obj.name: obj for obj in self._params}

    @property
    def names(self):
        return tuple([obj.name for obj in self._params])

    @property
    def labels(self):
        return tuple([obj.tex for obj in self._params])

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
    def tags(self):
        return tuple([obj.tag for obj in self._params])

    @property
    def params(self):
        return self._params

    def to_dict(self):
        return {obj.name: obj.value for obj in self._params}

    def from_tag(self, tag, values=False, index=False, invert=False):
        if values and index: assert 0
        tag = np.atleast_1d(tag)
        if not invert:
            ps = [(idx, obj) for idx, obj in enumerate(self._params)
                  if obj.tag in tag]
        else:
            ps = [(idx, obj) for idx, obj in enumerate(self._params)
                  if obj.tag not in tag]
        if values:
            return tuple([io[1].value for io in ps])
        elif index:
            return tuple([io[0] for io in ps])
        else:
            return ParamSet([io[1] for io in ps])


class SortingHelpFormatter(argparse.HelpFormatter):
    """Sort argparse help options alphabetically."""
    def add_arguments(self, actions):
        actions = sorted(actions, key=attrgetter('option_strings'))
        super(SortingHelpFormatter, self).add_arguments(actions)


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


def thread_type(t):
    if t.lower() == 'max':
        return multiprocessing.cpu_count()
    else:
        return int(t)


def thread_factors(t):
    for x in reversed(range(int(np.ceil(np.sqrt(t)))+1)):
        if t%x == 0:
            return (x, int(t/x))

