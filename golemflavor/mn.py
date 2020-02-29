# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : April 19, 2018

"""
Useful functions to use MultiNest for the BSM flavor ratio analysis
"""

from __future__ import absolute_import, division, print_function

from functools import partial

import numpy as np

from pymultinest import analyse, run

from golemflavor import llh as llh_utils
from golemflavor.misc import gen_identifier, make_dir, parse_bool, solve_ratio


def CubePrior(cube, ndim, n_params):
    pass


def lnProb(cube, ndim, n_params, mn_paramset, llh_paramset, asimov_paramset,
           args):
    if ndim != len(mn_paramset):
        raise AssertionError(
            'Length of MultiNest scan paramset is not the same as the input '
            'params\ncube={0}\nmn_paramset]{1}'.format(cube, mn_paramset)
        )
    pranges = mn_paramset.ranges
    for i in range(ndim):
        mn_paramset[i].value = (pranges[i][1]-pranges[i][0])*cube[i] + pranges[i][0]
    for pm in mn_paramset.names:
        llh_paramset[pm].value = mn_paramset[pm].value
    theta = llh_paramset.values
    llh = llh_utils.ln_prob(
        theta=theta,
        args=args,
        asimov_paramset=asimov_paramset,
        llh_paramset=llh_paramset,
    )
    return llh


def mn_argparse(parser):
    parser.add_argument(
        '--mn-live-points', type=int, default=3000,
        help='Number of live points for MultiNest evaluations'
    )
    parser.add_argument(
        '--mn-tolerance', type=float, default=0.01,
        help='Tolerance for MultiNest'
    )
    parser.add_argument(
        '--mn-efficiency', type=float, default=0.3,
        help='Sampling efficiency for MultiNest'
    )
    parser.add_argument(
        '--mn-output', type=str, default='./mnrun/',
        help='Folder to store MultiNest evaluations'
    )
    parser.add_argument(
        '--run-mn', type=parse_bool, default='True',
        help='Run MultiNest'
    )


def mn_evidence(mn_paramset, llh_paramset, asimov_paramset, args, prefix='mn'):
    """Run the MultiNest algorithm to calculate the evidence."""
    n_params = len(mn_paramset)

    for n in mn_paramset.names:
        llh_paramset[n].value = mn_paramset[n].value

    lnProbEval = partial(
        lnProb,
        mn_paramset     = mn_paramset,
        llh_paramset    = llh_paramset,
        asimov_paramset = asimov_paramset,
        args            = args,
    )

    if args.run_mn:
        make_dir(prefix)
        print('Running evidence calculation for {0}'.format(prefix))
        run(
            LogLikelihood              = lnProbEval,
            Prior                      = CubePrior,
            n_dims                     = n_params,
            n_live_points              = args.mn_live_points,
            evidence_tolerance         = args.mn_tolerance,
            sampling_efficiency        = args.mn_efficiency,
            outputfiles_basename       = prefix,
            importance_nested_sampling = True,
            # resume                     = False,
            resume                     = True,
            verbose                    = True
        )

    analyser = analyse.Analyzer(
        outputfiles_basename=prefix, n_params=n_params
    )
    return analyser.get_stats()['global evidence']
