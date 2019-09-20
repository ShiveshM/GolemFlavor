# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : March 17, 2018

"""
Useful GolemFit wrappers for the BSM flavour ratio analysis
"""

from __future__ import absolute_import, division

from functools import partial

import numpy as np

try:
    import GolemFitPy as gf
except:
    print 'Running without GolemFit'
    pass

from utils.enums import DataType, Likelihood, SteeringCateg
from utils.misc import enum_parse, parse_bool, thread_factors
from utils.param import ParamSet


FITTER = None


def fit_flags(llh_paramset):
    default_flags = {
        # False means it's not fixed in minimization
        'astroFlavorAngle1'         : True,
        'astroFlavorAngle2'         : True,
        'convNorm'                  : True,
        'promptNorm'                : True,
        'muonNorm'                  : True,
        'astroNorm'                 : True,
        'astroParticleBalance'      : True,
        # 'astroDeltaGamma'           : True,
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
    gf_nuisance = []
    for param in llh_paramset:
        if param.name in default_flags:
            print 'Setting param {0:<15} to float in the ' \
                'minimisation'.format(param.name)
            flags.__setattr__(param.name, False)
            gf_nuisance.append(param)
    return flags, ParamSet(gf_nuisance)


def steering_params(args):
    steering_categ = args.ast
    params = gf.SteeringParams(gf.sampleTag.MagicTau)
    params.quiet = False
    if args.debug:
        params.fastmode = False
    else:
        params.fastmode = True
    params.simToLoad= steering_categ.name.lower()
    params.evalThreads = args.threads

    if hasattr(args, 'binning'):
        params.minFitEnergy = args.binning[0]  # GeV
        params.maxFitEnergy = args.binning[-1] # GeV
    else:
        params.minFitEnergy = 6E4 # GeV
        params.maxFitEnergy = 1E7 # GeV
    params.load_data_from_text_file = False
    params.do_HESE_reshuffle=False
    params.use_legacy_selfveto_calculation = False

    return params


def setup_asimov(params):
    print 'Injecting the model', params
    asimov_params = gf.FitParameters(gf.sampleTag.MagicTau)
    for parm in params:
        asimov_params.__setattr__(parm.name, float(parm.value))
    FITTER.SetupAsimov(asimov_params)


def setup_realisation(params, seed):
    print 'Injecting the model', params
    asimov_params = gf.FitParameters(gf.sampleTag.MagicTau)
    for parm in params:
        asimov_params.__setattr__(parm.name, float(parm.value))
    FITTER.Swallow(FITTER.SpitRealization(asimov_params, seed))


def setup_fitter(args, asimov_paramset):
    global FITTER
    datapaths = gf.DataPaths()
    sparams = steering_params(args)
    npp = gf.NewPhysicsParams()
    FITTER = gf.GolemFit(datapaths, sparams, npp)
    if args.data is DataType.ASIMOV:
        setup_asimov(FITTER, asimov_paramset)
    elif args.data is DataType.REALISATION:
        seed = args.seed if args.seed is not None else 1
        setup_realisation(FITTER, asimov_paramset, seed)
    elif args.data is DataType.REAL:
        print 'Using MagicTau DATA'


def get_llh(params):
    fitparams = gf.FitParameters(gf.sampleTag.MagicTau)
    for parm in params:
        fitparams.__setattr__(parm.name, float(parm.value))
    llh = -FITTER.EvalLLH(fitparams)
    return llh


def get_llh_freq(params):
    print 'setting to {0}'.format(params)
    fitparams = gf.FitParameters(gf.sampleTag.MagicTau)
    for parm in params:
        fitparams.__setattr__(parm.name, float(parm.value))
    FITTER.SetFitParametersSeed(fitparams)
    llh = -FITTER.MinLLH().likelihood
    return llh


def data_distributions():
    hdat = FITTER.GetDataDistribution()
    binedges = np.asarray(
        [FITTER.GetZenithBinsData(), FITTER.GetEnergyBinsData()]
    )
    return hdat, binedges


def gf_argparse(parser):
    parser.add_argument(
        '--debug', default='False', type=parse_bool, help='Run without fastmode'
    )
    parser.add_argument(
        '--data', default='asimov', type=partial(enum_parse, c=DataType),
        choices=DataType, help='select datatype'
    )
    parser.add_argument(
        '--ast', default='p2_0', type=partial(enum_parse, c=SteeringCateg),
        choices=SteeringCateg,
        help='use asimov/fake dataset with specific steering'
    )
