# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : March 17, 2018

"""
Useful GolemFit wrappers for the BSM flavour ratio analysis
"""

from __future__ import absolute_import, division

import socket
from functools import partial

import GolemFitPy as gf

from utils.enums import DataType, SteeringCateg
from utils.misc import enum_parse, thread_factors


def fit_flags(llh_paramset):
    default_flags = {
        # False means it's not fixed in minimization
        'astroFlavorAngle1'         : True,
        'astroFlavorAngle2'         : True,
        'astroENorm'                : True,
        'astroMuNorm'               : True,
        'astroTauNorm'              : True,
        'convNorm'                  : True,
        'promptNorm'                : True,
        'muonNorm'                  : True,
        'astroNorm'                 : True,
        'astroParticleBalance'      : True,
        'astroDeltaGamma'           : True,
        'cutoffEnergy'              : True,
        'CRDeltaGamma'              : True,
        'piKRatio'                  : True,
        'NeutrinoAntineutrinoRatio' : True,
        'darkNorm'                  : True,
        'domEfficiency'             : True,
        'holeiceForward'            : True,
        'anisotropyScale'           : True,
        'astroNormSec'              : True,
        'astroDeltaGammaSec'        : True
    }
    flags = gf.FitParametersFlag()
    for param in llh_paramset:
        flags.__setattr__(param.name, False)
    return flags


def steering_params(args):
    steering_categ = args.ast
    params = gf.SteeringParams()
    params.quiet = False
    params.fastmode = True
    params.simToLoad= steering_categ.name.lower()
    params.spline_dom_efficiency = False
    params.spline_hole_ice = False
    params.spline_anisotrophy = False
    params.evalThreads = args.threads
    # params.evalThreads = thread_factors(args.threads)[1]
    params.diffuse_fit_type = gf.DiffuseFitType.SinglePowerLaw
    return params


def set_up_as(fitter, params):
    print 'Injecting the model', params
    asimov_params = gf.FitParameters(gf.sampleTag.HESE)
    for parm in params:
        asimov_params.__setattr__(parm.name, parm.value)
    fitter.SetupAsimov(asimov_params)


def setup_fitter(args, asimov_paramset):
    datapaths = gf.DataPaths()
    sparams = steering_params(args)
    npp = gf.NewPhysicsParams()
    fitter = gf.GolemFit(datapaths, sparams, npp)
    set_up_as(fitter, asimov_paramset)
    return fitter


def get_llh(fitter, params):
    fitparams = gf.FitParameters(gf.sampleTag.HESE)
    for parm in params:
        fitparams.__setattr__(parm.name, parm.value)
    llh = -fitter.EvalLLH(fitparams)
    return llh


def get_llh_freq(fitter, params):
    print 'setting to {0}'.format(params)
    fitparams = gf.FitParameters(gf.sampleTag.HESE)
    for parm in params:
        fitparams.__setattr__(parm.name, parm.value)
    fitter.SetFitParametersSeed(fitparams)
    return -fitter.MinLLH().likelihood


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

