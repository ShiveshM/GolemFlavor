# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : April 19, 2018

"""
Param class and functions for the BSM flavour ratio analysis
"""

from __future__ import absolute_import, division

import sys

from collections import Sequence
from copy import deepcopy

import numpy as np

from utils.plot import get_units
from utils.fr import fr_to_angles
from utils.enums import DataType, Likelihood, MixingScenario, ParamTag, PriorsCateg


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


def get_paramsets(args, nuisance_paramset):
    """Make the paramsets for generating the Asmimov MC sample and also running
    the MCMC.
    """
    asimov_paramset = []
    llh_paramset = []

    llh_paramset.extend(
        [x for x in nuisance_paramset.from_tag(ParamTag.SM_ANGLES)]
    )
    if args.likelihood in [Likelihood.GOLEMFIT, Likelihood.GF_FREQ]:
        gf_nuisance = [x for x in nuisance_paramset.from_tag(ParamTag.NUISANCE)]
        asimov_paramset.extend(gf_nuisance)
        llh_paramset.extend(gf_nuisance)
    for parm in llh_paramset:
        parm.value = args.__getattribute__(parm.name)
    tag = ParamTag.BESTFIT
    try:
        flavour_angles = fr_to_angles(args.measured_ratio)
    except:
        flavour_angles = fr_to_angles(args.injected_ratio)
    asimov_paramset.extend([
        Param(name='astroFlavorAngle1', value=flavour_angles[0], ranges=[0., 1.], std=0.2, tag=tag),
        Param(name='astroFlavorAngle2', value=flavour_angles[1], ranges=[-1., 1.], std=0.2, tag=tag),
    ])
    asimov_paramset = ParamSet(asimov_paramset)

    if hasattr(args, 'fix_source_ratio'):
        if args.fix_mixing is MixingScenario.NONE and not args.fix_mixing_almost:
            tag = ParamTag.MMANGLES
            llh_paramset.extend([
                Param(name='np_s_12^2', value=0.5, ranges=[0., 1.], std=0.2, tex=r'\tilde{s}_{12}^2', tag=tag),
                Param(name='np_c_13^4', value=0.5, ranges=[0., 1.], std=0.2, tex=r'\tilde{c}_{13}^4', tag=tag),
                Param(name='np_s_23^2', value=0.5, ranges=[0., 1.], std=0.2, tex=r'\tilde{s}_{23}^2', tag=tag),
                Param(name='np_dcp', value=np.pi, ranges=[0., 2*np.pi], std=0.2, tex=r'\tilde{\delta_{CP}}', tag=tag)
            ])
        if args.fix_mixing_almost:
            tag = ParamTag.MMANGLES
            llh_paramset.extend([
                Param(name='np_s_23^2', value=0.5, ranges=[0., 1.], std=0.2, tex=r'\tilde{s}_{23}^4', tag=tag)
            ])
        if not args.fix_scale:
            tag = ParamTag.SCALE
            if hasattr(args, 'dimension'):
                llh_paramset.append(
                    Param(name='logLam', value=np.log10(args.scale), ranges=np.log10(args.scale_region), std=3,
                          tex=r'{\rm log}_{10}\left (\Lambda^{-1}'+get_units(args.dimension)+r'\right )', tag=tag)
                )
            elif hasattr(args, 'dimensions'):
                llh_paramset.append(
                    Param(name='logLam', value=np.log10(args.scale), ranges=np.log10(args.scale_region), std=3,
                          tex=r'{\rm log}_{10}\left (\Lambda^{-1} / GeV^{-d+4}\right )', tag=tag)
                )
        if not args.fix_source_ratio:
            tag = ParamTag.SRCANGLES
            llh_paramset.extend([
                Param(name='s_phi4', value=0.5, ranges=[0., 1.],  std=0.2, tex=r'sin^4(\phi)', tag=tag),
                Param(name='c_2psi', value=0.5, ranges=[-1., 1.], std=0.2, tex=r'cos(2\psi)', tag=tag)
            ])
    llh_paramset = ParamSet(llh_paramset)
    return asimov_paramset, llh_paramset
