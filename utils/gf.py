# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : March 17, 2018

"""
Useful GolemFit wrappers for the BSM flavour ratio analysis
"""

from __future__ import absolute_import, division

import argparse
import socket
from functools import partial

import GolemFitPy as gf

from utils.enums import DataType, SteeringCateg
from utils.misc import enum_parse, thread_factors


def steering_params(args):
    steering_categ = args.ast
    params = gf.SteeringParams()
    if 'cobalt0' in socket.gethostname().split('.')[0]:
        params.quiet = False
    else:
        params.quiet = True
    # TODO(shivesh): figure out right number for missid
    params.track_to_shower_missid = 0.3
    params.fastmode = True
    # params.fastmode = False
    # params.readCompact = True
    params.readCompact = False
    params.simToLoad= steering_categ.name.lower()
    params.evalThreads = args.threads
    # params.evalThreads = thread_factors(args.threads)[1]
    params.spline_hole_ice = True
    params.spline_dom_efficiency = True
    return params


def set_up_as(fitter, params):
    print 'Injecting the model', params
    asimov_params = gf.FitParameters(gf.sampleTag.HESE)
    for parm in params:
        asimov_params.__setattr__(parm.name, parm.value)
    fitter.SetupAsimov(asimov_params)


def get_llh(fitter, params):
    fitparams = gf.FitParameters(gf.sampleTag.HESE)
    for parm in params:
        fitparams.__setattr__(parm.name, parm.value)
    return fitter.EvalLLH(fitparams)


def data_distributions(fitter):
    hdat = fitter.GetDataDistribution()
    binedges = np.asarray([fitter.GetZenithBinsData(), fitter.GetEnergyBinsData()])
    return hdat, binedges


def gf_argparse(parser):
    parser.add_argument(
        '--data', default='real', type=partial(enum_parse, c=DataType),
        choices=DataType, help='select datatype'
    )
    parser.add_argument(
        '--ast', default='p2_0', type=partial(enum_parse, c=SteeringCateg),
        choices=SteeringCateg,
        help='use asimov/fake dataset with specific steering'
    )

