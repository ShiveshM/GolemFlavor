#! /usr/bin/env python
# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : April 28, 2018

"""
HESE BSM flavour ratio analysis plotting script
"""

from __future__ import absolute_import, division

import os
import argparse
from functools import partial
from copy import deepcopy

import numpy as np
import numpy.ma as ma

from utils import fr as fr_utils
from utils import llh as llh_utils
from utils import misc as misc_utils
from utils import plot as plot_utils
from utils.enums import EnergyDependance, Likelihood, MixingScenario, ParamTag
from utils.enums import PriorsCateg, SensitivityCateg, StatCateg
from utils.param import Param, ParamSet


def process_args(args):
    """Process the input args."""
    if args.data is not DataType.REAL:
        args.injected_ratio = fr_utils.normalise_fr(args.injected_ratio)

    if len(args.source_ratios) % 3 != 0:
        raise ValueError(
            'Invalid source ratios input {0}'.format(args.source_ratios)
        )

    srs = [args.source_ratios[3*x:3*x+3]
           for x in range(int(len(args.source_ratios)/3))]
    args.source_ratios = map(fr_utils.normalise_fr, srs)

    args.dimensions = np.sort(args.dimensions)


def parse_args(args=None):
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="HESE BSM flavour ratio analysis plotting script",
        formatter_class=misc_utils.SortingHelpFormatter,
    )
    parser.add_argument(
        '--datadir', type=str,
        default='/data/user/smandalia/flavour_ratio/data/sensitivity',
        help='Path to data directory'
    )
    parser.add_argument(
        '--segments', type=int, default=10,
        help='Number of new physics scales to evaluate'
    )
    parser.add_argument(
        '--split-jobs', type=misc_utils.parse_bool, default='True',
        help='Did the jobs get split'
    )
    parser.add_argument(
        '--dimensions', type=int, nargs='*', default=[3, 6],
        help='Set the new physics dimensions to consider'
    )
    parser.add_argument(
        '--source-ratios', type=int, nargs='*', default=[1, 2, 0],
        help='Set the source flavour ratios for the case when you want to fix it'
    )
    parser.add_argument(
        '--texture', type=partial(enum_parse, c=Texture),
        default='none', choices=Texture, help='Set the BSM mixing texture'
    )
    parser.add_argument(
        '--data', default='asimov', type=partial(enum_parse, c=DataType),
        choices=DataType, help='select datatype'
    )
    parser.add_argument(
        '--plot-x', type=misc_utils.parse_bool, default='True',
        help='Make sensitivity plot x vs limit'
    )
    parser.add_argument(
        '--plot-table', type=misc_utils.parse_bool, default='True',
        help='Make sensitivity table plot'
    )
    parser.add_argument(
        '--plot-statistic', type=misc_utils.parse_bool, default='False',
        help='Plot MultiNest evidence or LLH value'
    )
    llh_utils.likelihood_argparse(parser)
    if args is None: return parser.parse_args()
    else: return parser.parse_args(args.split())


def main():
    args = parse_args()
    process_args(args)
    misc_utils.print_args(args)

    dims = len(args.dimensions)
    srcs = len(args.source_ratios)
    if args.texture is Texture.NONE:
        textures = [Texture.OEU, Texture.OET, Texture.OUT]
    else:
        textures = [args.texture]
    texs = len(textures)

    prefix = ''
    # prefix = 'noprior'

    # Initialise data structure.
    statistic_arr = np.full((dims, srcs, texs, args.segments, 2), np.nan)

    print 'Loading data'
    for idim, dim in enumerate(args.dimensions):
        argsc = deepcopy(args)
        argsc.dimension = dim

        # Array of scales to scan over.
        boundaries = fr_utils.SCALE_BOUNDARIES[argsc.dimension]
        eval_scales = np.linspace(
            boundaries[0], boundaries[1], args.segments-1
        )
        eval_scales = np.concatenate([[-100.], eval_scales])

        for isrc, src in enumerate(args.source_ratios):
            argsc.source_ratio = src
            for itex, texture in enumerate(textures):
                argc.texture = texture

                base_infile = args.datadir + '/{0}/{1}/{2}/fr_stat'.format(
                    *map(misc_utils.parse_enum, [args.stat_method, args.data]),
                    prefix
                ) + misc_utils.gen_identifier(argsc)

                print '== {0:<25} = {1}'.format('base_infile', base_infile)

                if args.split_jobs:
                    for idx_sc, scale in enumerate(eval_scales):
                        infile = base_infile + '_scale_{0:.0E}'.format(
                            np.power(10, scale)
                        )
                        try:
                            print 'Loading from {0}'.format(infile+'.npy')
                            statistic_arr[idim][isrc][itex][idx_sc] = \
                                np.load(infile+'.npy')[:,0]
                        except:
                            print 'Unable to load file {0}'.format(
                                infile+'.npy'
                            )
                            continue
                else:
                    print 'Loading from {0}'.format(base_infile+'.npy')
                    try:
                        statistic_arr[idim][isrc][itex] = \
                            np.load(base_infile+'.npy')
                    except:
                        print 'Unable to load file {0}'.format(
                            base_infile+'.npy'
                        )
                        continue

    data = ma.masked_invalid(statistic_arr)

    print 'data', data
    if args.plot_statistic:
        print 'Plotting statistic'

        for idim, dim in enumerate(args.dimensions):
            argsc = deepcopy(args)
            argsc.dimension = dim

            # Array of scales to scan over.
            boundaries = fr_utils.SCALE_BOUNDARIES[argsc.dimension]
            eval_scales = np.linspace(
                boundaries[0], boundaries[1], args.segments-1
            )
            eval_scales = np.concatenate([[-100.], eval_scales])

            for isrc, src in enumerate(args.source_ratios):
                argsc.source_ratio = src
                for itex, texture in enumerate(textures):
                    argc.texture = texture

                    base_infile = args.datadir + '/{0}/{1}/{2}/fr_stat'.format(
                        *map(misc_utils.parse_enum, [args.stat_method, args.data]),
                        prefix
                    ) + misc_utils.gen_identifier(argsc)
                    basename = os.path.dirname(base_infile)
                    outfile = basename[:5]+basename[5:].replace('data', 'plots')
                    outfile += '/' + os.path.basename(base_infile)

                    label = plot_utils.texture_label(texture)[:-1]+r'=\pi/4$'
                    plot_utils.plot_statistic(
                        data        = data[idim][isrc][itex],
                        outfile     = outfile,
                        outformat   = ['png'],
                        args        = argsc,
                        scale_param = scale,
                        label       = label
                    )

    basename = args.datadir[:5]+args.datadir[5:].replace('data', 'plots')
    baseoutfile = basename + '/{0}/{1}/{2}/'.format(
        *map(misc_utils.parse_enum, [args.stat_method, args.data]), prefix
    )

    if args.plot_x:
        for idim, dim in enumerate(args.dimensions):
            argsc = deepcopy(args)
            argsc.dimension = dim
            plot_utils.plot_x(
                data      = data,
                outfile   = baseoutfile + '/hese_x',
                outformat = ['png', 'pdf'],
                args      = argsc,
            )

    if args.plot_table:
        plot_utils.plot_table_sens(
            data      = data,
            outfile   = baseoutfile + '/hese_table',
            outformat = ['png', 'pdf'],
            args      = args,
        )


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
