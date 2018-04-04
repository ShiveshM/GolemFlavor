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

from utils.enums import *
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
    params.do_HESE_reshuffle = False
    params.sampleToLoad = gf.sampleTag.HESE
    params.simToLoad= steering_categ.name.lower()
    params.evalThreads = args.threads
    # params.evalThreads = thread_factors(args.threads)[1]
    params.baseline_astro_spectral_index = -2.
    params.spline_hole_ice = True
    params.spline_dom_efficiency = True
    if steering_categ == SteeringCateg.LONGLIFE:
        params.years = [999]
        params.simToLoad= 'p2_0'
    elif steering_categ == SteeringCateg.DPL:
        params.diffuse_fit_type = gf.DiffuseFitType.DoublePowerLaw
        params.simToLoad= 'p2_0'
    return params


def set_up_as(fitter, params):
    print 'Injecting the model', params
    asimov_params = gf.FitParameters()
    for parm in params:
        asimov_params.__setattr__(parm.name, parm.value)
    fitter.SetupAsimov(asimov_params)


def get_llh(fitter, params):
    fitparams = gf.FitParameters()
    for parm in params:
        fitparams.__setattr__(parm.name, parm.value)
    return fitter.EvalLLH(fitparams)


def data_distributions(fitter):
    hdat = fitter.GetDataDistribution()
    binedges = np.asarray([fitter.GetZenithBinsData(), fitter.GetEnergyBinsData()])
    return hdat, binedges


def fit_flags(fitflag_categ):
    flags = gf.FitParametersFlag()
    if fitflag_categ is FitFlagCateg.xs:
        # False means it's not fixed in minimization
        flags.NeutrinoAntineutrinoRatio = False
    return flags


def fit_params(fit_categ):
    params = gf.FitParameters()
    params.astroNorm = 7.5
    params.astroDeltaGamma = 0.9
    if fit_categ is FitCateg.hesespl:
        params.astroNormSec = 0
    elif fit_categ is FitCateg.hesedpl:
        params.astroNormSec=7.0
    elif fit_categ is FitCateg.zpspl:
        # zero prompt, single powerlaw
        params.promptNorm = 0
        params.astroNormSec = 0
    elif fit_categ is FitCateg.zpdpl:
        # zero prompt, double powerlaw
        params.promptNorm = 0
        params.astroNormSec=7.0
    elif fit_categ is FitCateg.nunubar2:
        params.NeutrinoAntineutrinoRatio = 2


def fit_priors(fitpriors_categ):
    priors = gf.Priors()
    if fitpriors_categ == FitPriorsCateg.xs:
        priors.promptNormCenter = 1
        priors.promptNormWidth = 3
        priors.astroDeltaGammaCenter = 0
        priors.astroDeltaGammaWidth = 1
    return priors


def gen_steering_params(steering_categ, quiet=False):
    params = gf.SteeringParams()
    if quiet: params.quiet = True
    params.fastmode = False
    params.do_HESE_reshuffle = False
    params.numc_tag = steering_categ.name
    params.baseline_astro_spectral_index = -2.
    if steering_categ is SteeringCateg.LONGLIFE:
        params.years = [999]
        params.numc_tag = 'std_half1'
    if steering_categ is SteeringCateg.DPL:
        params.diffuse_fit_type = gf.DiffuseFitType.DoublePowerLaw
        params.numc_tag = 'std_half1'
    return params


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
    parser.add_argument(
        '--aft', default='hesespl', type=partial(enum_parse, c=FitCateg),
        choices=FitCateg,
        help='use asimov/fake dataset with specific Fit'
    )
    parser.add_argument(
        '--axs', default='nom', type=partial(enum_parse, c=XSCateg),
        choices=XSCateg,
        help='use asimov/fake dataset with xs scaling'
    )
    parser.add_argument(
        '--priors', default='uniform', type=partial(enum_parse, c=Priors),
        choices=Priors, help='Bayesian priors'
    )

