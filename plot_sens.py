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
from utils import plot as plot_utils
from utils.enums import DataType, Texture
from utils.misc import enum_parse, parse_bool, parse_enum, print_args
from utils.misc import gen_identifier, SortingHelpFormatter
from utils.param import Param, ParamSet


MASK_X = (0.3, 0.7)


def process_args(args):
    """Process the input args."""
    if args.data is not DataType.REAL:
        args.injected_ratio = fr_utils.normalise_fr(args.injected_ratio)

    # Anon points
    anon = []
    if args.dimensions[0] == 3:
        anon.append([0.825, 0.845])
        anon.append([0.865, 0.875])
        anon.append([0.875, 0.885])
        anon.append([0.905, 0.915])
        anon.append([0.925, 0.935])
    if args.dimensions[0] == 4:
        anon.append([0.165, 0.175])
        anon.append([0.805, 0.825])
        anon.append([0.835, 0.845])
        anon.append([0.855, 0.885])
        anon.append([0.965, 0.975])
    if args.dimensions[0] == 5:
        anon.append([0.895, 0.905])
        anon.append([0.955, 0.965])
    if args.dimensions[0] == 6:
        anon.append([0.115, 0.125])
        anon.append([0.855, 0.865])
    if args.dimensions[0] == 7:
        # anon.append([0.815, 0.835])
        anon.append([0.875, 0.885])
    if args.dimensions[0] == 8:
        anon.append([0.915, 0.935])
        anon.append([0.875, 0.895])
        anon.append([0.845, 0.855])

    if args.source_ratios is not None:
        if args.x_segments is not None:
            raise ValueError('Cannot do both --source-ratios and --x-segments')
        if len(args.source_ratios) % 3 != 0:
            raise ValueError(
                'Invalid source ratios input {0}'.format(args.source_ratios)
            )

        srs = [args.source_ratios[3*x:3*x+3]
               for x in range(int(len(args.source_ratios)/3))]
        args.source_ratios = map(fr_utils.normalise_fr, srs)
    elif args.x_segments is not None:
        x_array = np.linspace(0, 1, args.x_segments)
        sources = []
        for x in x_array:
            if x >= MASK_X[0] and x <= MASK_X[1]: continue
            skip = False
            for a in anon:
                if (a[1] > x) & (x > a[0]):
                    print 'Skipping src', x
                    skip = True
                    break
            if skip: continue
            sources.append([x, 1-x, 0])
        args.source_ratios = sources
    else:
        raise ValueError('Must supply either --source-ratios or --x-segments')

    args.dimensions = np.sort(args.dimensions)


def parse_args(args=None):
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="HESE BSM flavour ratio analysis plotting script",
        formatter_class=SortingHelpFormatter,
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
        '--split-jobs', type=parse_bool, default='True',
        help='Did the jobs get split'
    )
    parser.add_argument(
        '--dimensions', type=int, nargs='*', default=[3, 6],
        help='Set the new physics dimensions to consider'
    )
    parser.add_argument(
        '--source-ratios', type=int, nargs='*', default=None,
        required=False, help='Set the source flavour ratios'
    )
    parser.add_argument(
        '--x-segments', type=int, default=None,
        required=False, help='Number of segments in x'
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
        '--plot-x', type=parse_bool, default='False',
        help='Make sensitivity plot x vs limit'
    )
    parser.add_argument(
        '--plot-table', type=parse_bool, default='False',
        help='Make sensitivity table plot'
    )
    parser.add_argument(
        '--plot-statistic', type=parse_bool, default='False',
        help='Plot MultiNest evidence or LLH value'
    )
    llh_utils.llh_argparse(parser)
    if args is None: return parser.parse_args()
    else: return parser.parse_args(args.split())


def main():
    args = parse_args()
    process_args(args)
    print_args(args)

    dims = len(args.dimensions)
    srcs = len(args.source_ratios)
    if args.texture is Texture.NONE:
        textures = [Texture.OEU, Texture.OET, Texture.OUT]
        if args.plot_table:
            textures = [Texture.OET, Texture.OUT]
    else:
        textures = [args.texture]
    texs = len(textures)

    prefix = ''
    # prefix = 'noprior'

    # Initialise data structure.
    statistic_arr = np.full((dims, srcs, texs, args.segments, 2), np.nan)

    print 'Loading data'
    argsc = deepcopy(args)
    for idim, dim in enumerate(args.dimensions):
        argsc.dimension = dim

        datadir = args.datadir + '/DIM{0}'.format(dim)
        # Array of scales to scan over.
        boundaries = fr_utils.SCALE_BOUNDARIES[argsc.dimension]
        eval_scales = np.linspace(
            boundaries[0], boundaries[1], args.segments-1
        )
        eval_scales = np.concatenate([[-100.], eval_scales])

        for isrc, src in enumerate(args.source_ratios):
            argsc.source_ratio = src

            for itex, texture in enumerate(textures):
                argsc.texture = texture

                base_infile = datadir + '/{0}/{1}/'.format(
                    *map(parse_enum, [args.stat_method, args.data])
                ) + r'{0}/fr_stat'.format(prefix) + gen_identifier(argsc)

                print '== {0:<25} = {1}'.format('base_infile', base_infile)

                if args.split_jobs:
                    for idx_sc, scale in enumerate(eval_scales):
                        infile = base_infile + '_scale_{0:.0E}'.format(
                            np.power(10, scale)
                        )
                        try:
                            print 'Loading from {0}'.format(infile+'.npy')
                            statistic_arr[idim][isrc][itex][idx_sc] = \
                                np.load(infile+'.npy')[0]
                        except:
                            print 'Unable to load file {0}'.format(
                                infile+'.npy'
                            )
                            # raise
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
                    argsc.texture = texture

                    base_infile = args.datadir + '/DIM{0}/{1}/{2}/'.format(
                        dim, *map(parse_enum, [args.stat_method, args.data])
                    ) + r'{0}/fr_stat'.format(prefix) + gen_identifier(argsc)
                    basename = os.path.dirname(base_infile)
                    outfile = basename[:5]+basename[5:].replace('data', 'plots')
                    outfile += '/' + os.path.basename(base_infile)

                    label = r'$\text{Texture}=' + plot_utils.texture_label(texture)[1:]
                    plot_utils.plot_statistic(
                        data        = data[idim][isrc][itex],
                        outfile     = outfile,
                        outformat   = ['png'],
                        args        = argsc,
                        scale_param = scale,
                        label       = label
                    )

    basename = args.datadir[:5]+args.datadir[5:].replace('data', 'plots')
    baseoutfile = basename + '/{0}/{1}/'.format(
        *map(parse_enum, [args.stat_method, args.data])
    ) + r'{0}'.format(prefix)

    argsc = deepcopy(args)
    if args.plot_x:
        for idim, dim in enumerate(args.dimensions):
            print '|||| DIM = {0}'.format(dim)
            argsc.dimension = dim
            plot_utils.plot_x(
                data      = data[idim],
                outfile   = baseoutfile + '/hese_x_DIM{0}'.format(dim),
                outformat = ['png', 'pdf'],
                args      = argsc,
                normalise = True
            )

    if args.plot_table:
        plot_utils.plot_table_sens(
            data      = data,
            outfile   = baseoutfile + '/hese_table',
            outformat = ['png', 'pdf'],
            args      = args,
            show_lvatmo = True
        )


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
