# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : March 17, 2018

"""
Useful GolemFit wrappers for the BSM flavour ratio analysis
"""

from __future__ import absolute_import, division

import argparse
from functools import partial

import GolemFitPy as gf

from utils.enums import *
from utils.misc import enum_parse


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
        '--ast', default='baseline', type=partial(enum_parse, c=SteeringCateg),
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

