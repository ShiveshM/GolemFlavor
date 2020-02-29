# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : April 19, 2018

"""
Param class and functions for the BSM flavor ratio analysis
"""

from __future__ import absolute_import, division

import sys

from collections import Sequence
from copy import deepcopy

import numpy as np

from golemflavor.fr import fr_to_angles
from golemflavor.enums import DataType, Likelihood, ParamTag, PriorsCateg


class Param(object):
    """Parameter class to store parameters."""
    def __init__(self, name, value, ranges, prior=None, seed=None, std=None,
                 tex=None, tag=None):
        self._prior = None
        self._seed = None
        self._ranges = None
        self._tex = None
        self._tag = None

        self.name = name
        self.value = value
        self.nominal_value = deepcopy(value)
        self.prior = prior
        self.ranges = ranges
        self.seed = seed
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
    def prior(self):
        return self._prior

    @prior.setter
    def prior(self, value):
        if value is None:
            self._prior = PriorsCateg.UNIFORM
        else:
            assert value in PriorsCateg
            self._prior = value

    @property
    def seed(self):
        if self._seed is None: return self.ranges
        return tuple(self._seed)

    @seed.setter
    def seed(self, values):
        if values is None: return
        self._seed = [val for val in values]

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
    """Container class for a set of parameters."""
    def __init__(self, *args):
        param_sequence = []
        for arg in args:
            try:
                param_sequence.extend(arg)
            except TypeError:
                param_sequence.append(arg)

        if len(param_sequence) != 0:
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
        return super(ParamSet, self).__getattribute__(attr)

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
    def nominal_values(self):
        return tuple([obj.nominal_value for obj in self._params])

    @property
    def seeds(self):
        return tuple([obj.seed for obj in self._params])

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

    def remove_params(self, params):
        rm_paramset = []
        for parm in self.params:
            if parm.name not in params.names:
                rm_paramset.append(parm)
        return ParamSet(rm_paramset)

    def extend(self, p):
        param_sequence = self.params
        if isinstance(p, Param):
            param_sequence.append(p)
        elif isinstance(p, ParamSet):
            param_sequence.extend(p.params)
        return ParamSet(param_sequence)
