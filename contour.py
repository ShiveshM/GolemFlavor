#! /usr/bin/env python
# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : November 26, 2018

"""
HESE flavour ratio contour
"""

from __future__ import absolute_import, division

import argparse
from functools import partial

import numpy as np

from utils import fr as fr_utils
from utils import gf as gf_utils
from utils import llh as llh_utils
from utils import misc as misc_utils
from utils import mn as mn_utils
from utils import plot as plot_utils
from utils.enums import str_enum
from utils.enums import DataType, Likelihood, ParamTag
from utils.param import Param, ParamSet, get_paramsets

from pymultinest import Analyzer, run


def define_nuisance():
    """Define the nuisance parameters."""
    nuisance = []
    tag = ParamTag.NUISANCE
    nuisance.extend([
        Param(name='convNorm',        value=1.,  seed=[0.5, 2. ], ranges=[0. , 50.], std=0.3,  tag=tag),
        Param(name='promptNorm',      value=0.,  seed=[0. , 6. ], ranges=[0. , 50.], std=0.05, tag=tag),
        Param(name='muonNorm',        value=1.,  seed=[0.1, 2. ], ranges=[0. , 50.], std=0.1,  tag=tag),
        Param(name='astroNorm',       value=6.9, seed=[0.1, 10.], ranges=[0. , 50.], std=0.1,  tag=tag),
        Param(name='astroDeltaGamma', value=2.5, seed=[2.4, 3. ], ranges=[-5., 5. ], std=0.1,  tag=tag)
    ])
    return ParamSet(nuisance)


def nuisance_argparse(parser):
    nuisance = define_nuisance()
    for parm in nuisance:
        parser.add_argument(
            '--'+parm.name, type=float, default=parm.value,
            help=parm.name+' to inject'
        )

def process_args(args):
    """Process the input args."""
    args.plot_angles = args.plot_chains
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
        description="BSM flavour ratio analysis",
        formatter_class=misc_utils.SortingHelpFormatter,
    )
    parser.add_argument(
        '--injected-ratio', type=float, nargs=3, default=[1, 1, 1],
        help='Set the central value for the injected flavour ratio at IceCube'
    )
    parser.add_argument(
        '--run-scan', type=misc_utils.parse_bool, default='True',
        help='Do the scan from scratch'
    )
    parser.add_argument(
        '--plot-chains', type=misc_utils.parse_bool, default='False',
        help='Plot the (joint) posteriors'
    )
    parser.add_argument(
        '--plot-triangle', type=misc_utils.parse_bool, default='False',
        help='Project the posterior contour on the flavour triangle'
    )
    parser.add_argument(
        '--seed', type=misc_utils.seed_parse, default='25',
        help='Set the random seed value'
    )
    parser.add_argument(
        '--threads', type=misc_utils.thread_type, default='1',
        help='Set the number of threads to use (int or "max")'
    )
    parser.add_argument(
        '--outfile', type=str, default='./untitled',
        help='Path to output results'
    )
    try:
        gf_utils.gf_argparse(parser)
    except: pass
    llh_utils.likelihood_argparse(parser)
    mn_utils.mn_argparse(parser)
    nuisance_argparse(parser)
    misc_utils.remove_option(parser, 'sigma_ratio')
    misc_utils.remove_option(parser, 'mn_output')
    if args is None: return parser.parse_args()
    else: return parser.parse_args(args.split())


def gen_identifier(args):
    f = '_{0}_{1}'.format(*map(str_enum, (args.likelihood, args.data)))
    if args.data is not DataType.REAL:
        ir1, ir2, ir3 = misc_utils.solve_ratio(args.injected_ratio)
        f += '_INJ_{0:03d}_{1:03d}_{2:03d}'.format(ir1, ir2, ir3)
    return f


def gen_figtext(args, asimov_paramset):
    f = ''
    if args.data is DataType.REAL:
        f += 'IceCube Preliminary\n'
    else:
        ir1, ir2, ir3 = misc_utils.solve_ratio(args.injected_ratio)
        f += 'Injected ratio = [{0}, {1}, {2}]\n'.format(ir1, ir2, ir3)
        for param in asimov_paramset:
            f += 'Injected {0} = {1:.3f}\n'.format(
                param.name, param.nominal_value
            )
    return f


def triangle_llh(theta, args, hypo_paramset, fitter):
    """Log likelihood function for a given theta."""
    if len(theta) != len(hypo_paramset):
        raise AssertionError(
            'Dimensions of scan is not the same as the input '
            'params\ntheta={0}\nparamset]{1}'.format(theta, hypo_paramset)
        )
    for idx, param in enumerate(hypo_paramset):
        param.value = theta[idx]

    if args.likelihood is Likelihood.GOLEMFIT:
        llh = gf_utils.get_llh(fitter, hypo_paramset)
    elif args.likelihood is Likelihood.GF_FREQ:
        llh = gf_utils.get_llh_freq(fitter, hypo_paramset)

    return llh


def ln_prob(theta, args, hypo_paramset, fitter):
    lp = llh_utils.lnprior(theta, paramset=hypo_paramset)
    if not np.isfinite(lp):
        return -np.inf
    return lp + triangle_llh(
        theta,
        args           = args,
        hypo_paramset  = hypo_paramset,
        fitter         = fitter
    )


def lnProb(cube, ndim, n_params, hypo_paramset, args, fitter):
    if ndim != len(hypo_paramset):
        raise AssertionError(
            'Length of MultiNest scan paramset is not the same as the input '
            'params\ncube={0}\nmn_paramset]{1}'.format(cube, hypo_paramset)
        )
    pranges = hypo_paramset.ranges
    for i in xrange(ndim):
        hypo_paramset[i].value = (pranges[i][1]-pranges[i][0])*cube[i] + pranges[i][0]
    theta = hypo_paramset.values
    llh = ln_prob(
        theta          = theta,
        args           = args,
        hypo_paramset  = hypo_paramset,
        fitter         = fitter
    )
    return llh


def main():
    args = parse_args()
    process_args(args)
    misc_utils.print_args(args)

    if args.seed is not None:
        np.random.seed(args.seed)

    asimov_paramset, hypo_paramset = get_paramsets(args, define_nuisance())
    hypo_paramset.extend(asimov_paramset.from_tag(ParamTag.BESTFIT))
    outfile = args.outfile + gen_identifier(args)
    print '== {0:<25} = {1}'.format('outfile', outfile)

    n_params = len(hypo_paramset)
    prefix = outfile + '_mn_'
    misc_utils.make_dir(prefix)

    print 'asimov_paramset', asimov_paramset
    print 'hypo_paramset', hypo_paramset

    if args.run_scan:
        fitter = gf_utils.setup_fitter(args, asimov_paramset)

        lnProbEval = partial(
            lnProb,
            hypo_paramset  = hypo_paramset,
            args           = args,
            fitter         = fitter
        )

        print 'Running evidence calculation for {0}'.format(prefix)
        run(
            LogLikelihood              = lnProbEval,
            Prior                      = mn_utils.CubePrior,
            n_dims                     = n_params,
            n_live_points              = args.mn_live_points,
            evidence_tolerance         = args.mn_tolerance,
            outputfiles_basename       = prefix,
            importance_nested_sampling = True,
            resume                     = False,
            verbose                    = True
        )

    # Analyze
    analyser = Analyzer(
        outputfiles_basename=prefix, n_params=n_params
    )
    print analyser

    pranges = hypo_paramset.ranges
    fig_text = gen_figtext(args, asimov_paramset)

    if args.plot_chains:
        chains = analyser.get_data()[:,2:]
        for x in chains:
            for i in xrange(len(x)):
                x[i] = (pranges[i][1]-pranges[i][0])*x[i] + pranges[i][0]
        plot_utils.chainer_plot(
            infile       = chains,
            outfile      = outfile[:5]+outfile[5:].replace('data', 'plots')+'_posterior',
            outformat    = ['pdf'],
            args         = args,
            llh_paramset = hypo_paramset,
            fig_text     = fig_text
        )

    if args.plot_triangle:
        llh = analyser.get_data()[:,1]

        chains = analyser.get_data()[:,2:]
        for x in chains:
            for i in xrange(len(x)):
                x[i] = (pranges[i][1]-pranges[i][0])*x[i] + pranges[i][0]
        flavour_angles = chains[:,-2:]
        flavour_ratios = np.array(
            map(fr_utils.angles_to_fr, flavour_angles), dtype=np.float
        )

        plot_utils.triangle_project(
            frs          = flavour_ratios,
            llh          = llh,
            outfile      = outfile[:5]+outfile[5:].replace('data', 'plots')+'_triangle',
            outformat    = ['png'],
            args         = args,
            llh_paramset = hypo_paramset,
            fig_text     = fig_text
        )

    print "DONE!"


main.__doc__ = __doc__


if __name__ == '__main__':
    main()
