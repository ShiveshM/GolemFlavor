#! /usr/bin/env python
# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : February 24, 2019

"""
HESE BSM Flavour Figure 2
"""

from __future__ import absolute_import, division

import argparse
from functools import partial

import numpy as np

from utils import fr as fr_utils
from utils import misc as misc_utils
from utils import plot as plot_utils
from utils.enums import str_enum
from utils.enums import Likelihood, ParamTag, PriorsCateg
from utils.param import Param, ParamSet

from matplotlib import pyplot as plt

# from pymultinest import Analyzer
import json


def define_nuisance():
    """Define the nuisance parameters."""
    nuisance = []
    tag = ParamTag.NUISANCE
    lg_prior = PriorsCateg.LIMITEDGAUSS
    nuisance.extend([
        Param(name='convNorm',        value=1.,  seed=[0.5, 2. ], ranges=[0.1, 10.], std=0.4, prior=lg_prior, tag=tag),
        Param(name='promptNorm',      value=0.,  seed=[0. , 6. ], ranges=[0. , 20.], std=2.4, prior=lg_prior, tag=tag),
        Param(name='muonNorm',        value=1.,  seed=[0.1, 2. ], ranges=[0. , 10.], std=0.1, tag=tag),
        Param(name='astroNorm',       value=6.9, seed=[0.,  5. ], ranges=[0. , 20.], std=1.5, tag=tag),
        Param(name='astroDeltaGamma', value=2.5, seed=[2.4, 3. ], ranges=[-5., 5. ], std=0.1, tag=tag)
    ])
    return ParamSet(nuisance)


def get_paramsets(args, nuisance_paramset):
    paramset = []
    if args.likelihood in [Likelihood.GOLEMFIT, Likelihood.GF_FREQ]:
        gf_nuisance = [x for x in nuisance_paramset.from_tag(ParamTag.NUISANCE)]
        paramset.extend(gf_nuisance)
    tag = ParamTag.BESTFIT
    paramset.extend([
        Param(name='astroFlavorAngle1', value=0, ranges=[0., 1.], std=0.2, tag=tag),
        Param(name='astroFlavorAngle2', value=0, ranges=[-1., 1.], std=0.2, tag=tag),
    ])
    paramset = ParamSet(paramset)
    return paramset


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
        description="HESE BSM Flavour Figure 2",
        formatter_class=misc_utils.SortingHelpFormatter,
    )
    parser.add_argument(
        '--likelihood', default='golemfit',
        type=partial(misc_utils.enum_parse, c=Likelihood),
        choices=Likelihood, help='likelihood contour'
    )
    parser.add_argument(
        '--contour-dir', type=str,
        help='Path to directory containing MultiNest runs'
    )
    parser.add_argument(
        '--outfile', type=str, default='./untitled',
        help='Output path'
    )
    if args is None: return parser.parse_args()
    else: return parser.parse_args(args.split())


def main():
    args = parse_args()
    process_args(args)
    misc_utils.print_args(args)

    paramset = get_paramsets(args, define_nuisance())
    n_params = len(paramset)
    print n_params

    # Data
    data_path = '/home/aschneider/programs/GOLEMSPACE/sources/GolemFit/scripts/diffuse/mcmcs/results/dpl_numu_prior_flavor_20190302-162221-a747f528-8aa6-4488-8c80-059572c099fe.json'
    with open(data_path) as f:
        d_json = json.load(f)

    names = d_json['func_args']
    chains = np.array(d_json['chain'])
    print 'names', names
    print 'chains.shape', chains.shape

    flavour_angles = chains[:,5:7]
    flavour_ratios = np.array(
        map(fr_utils.angles_to_fr, flavour_angles)
    )

    nbins = 25

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    tax = plot_utils.get_tax(ax, scale=nbins)

    plot_utils.flavour_contour(
        frs = flavour_ratios,
        ax = ax,
        nbins = nbins,
        coverage = 99,
        linewidth = 2,
        color = 'green'
    )

    plot_utils.flavour_contour(
        frs = flavour_ratios,
        ax = ax,
        nbins = nbins,
        coverage = 90,
        linewidth = 2,
        color = 'blue'
    )

    plot_utils.flavour_contour(
        frs = flavour_ratios,
        ax = ax,
        nbins = nbins,
        coverage = 68,
        linewidth = 2,
        color = 'red'
    )

    ax.legend()

    fig.savefig('test_austin.png', bbox_inches='tight', dpi=150)

    print "DONE!"


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
