# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : April 04, 2018

"""
Likelihood functions for the BSM flavour ratio analysis
"""

from __future__ import absolute_import, division

from copy import deepcopy
from functools import partial

import numpy as np
import scipy
from scipy.stats import multivariate_normal, truncnorm

from golemflavor import fr as fr_utils
from golemflavor import gf as gf_utils
from golemflavor.enums import Likelihood, ParamTag, PriorsCateg, StatCateg
from golemflavor.misc import enum_parse, gen_identifier, parse_bool


def GaussianBoundedRV(loc=0., sigma=1., lower=-np.inf, upper=np.inf):
    """Normalised gaussian bounded between lower and upper values"""
    low, up = (lower - loc) / sigma, (upper - loc) / sigma
    g = scipy.stats.truncnorm(loc=loc, scale=sigma, a=low, b=up)
    return g


def multi_gaussian(fr, fr_bf, sigma, offset=-320):
    """Multivariate gaussian likelihood."""
    cov_fr = np.identity(3) * sigma
    return np.log(multivariate_normal.pdf(fr, mean=fr_bf, cov=cov_fr)) + offset


def llh_argparse(parser):
    parser.add_argument(
        '--stat-method', default='bayesian',
        type=partial(enum_parse, c=StatCateg), choices=StatCateg,
        help='Statistical method to employ'
    )


def lnprior(theta, paramset):
    """Priors on theta."""
    if len(theta) != len(paramset):
        raise AssertionError(
            'Length of MCMC scan is not the same as the input '
            'params\ntheta={0}\nparamset={1}'.format(theta, paramset)
        )
    for idx, param in enumerate(paramset):
        param.value = theta[idx]
    ranges = paramset.ranges
    for value, range in zip(theta, ranges):
        if range[0] <= value <= range[1]:
            pass
        else: return -np.inf

    prior = 0
    for param in paramset:
        if param.prior is PriorsCateg.GAUSSIAN:
            prior += GaussianBoundedRV(
                loc=param.nominal_value, sigma=param.std
            ).logpdf(param.value)
        elif param.prior is PriorsCateg.LIMITEDGAUSS:
            prior += GaussianBoundedRV(
                loc=param.nominal_value, sigma=param.std,
                lower=param.ranges[0], upper=param.ranges[1]
            ).logpdf(param.value)
    return prior


def triangle_llh(theta, args, asimov_paramset, llh_paramset):
    """Log likelihood function for a given theta."""
    if len(theta) != len(llh_paramset):
        raise AssertionError(
            'Length of MCMC scan is not the same as the input '
            'params\ntheta={0}\nparamset]{1}'.format(theta, llh_paramset)
        )
    hypo_paramset = asimov_paramset
    for param in llh_paramset.from_tag(ParamTag.NUISANCE):
        hypo_paramset[param.name].value = param.value

    spectral_index = -hypo_paramset['astroDeltaGamma'].value
    # Assigning llh_paramset values from theta happens in this function.
    fr = fr_utils.flux_averaged_BSMu(theta, args, spectral_index, llh_paramset)

    flavour_angles = fr_utils.fr_to_angles(fr)
    # print 'flavour_angles', map(float, flavour_angles)
    for idx, param in enumerate(hypo_paramset.from_tag(ParamTag.BESTFIT)):
        param.value = flavour_angles[idx]

    if args.likelihood is Likelihood.GOLEMFIT:
        llh = gf_utils.get_llh(hypo_paramset)
    elif args.likelihood is Likelihood.GF_FREQ:
        llh = gf_utils.get_llh_freq(hypo_paramset)
    return llh


def ln_prob(theta, args, asimov_paramset, llh_paramset):
    dc_asimov_paramset = deepcopy(asimov_paramset)
    dc_llh_paramset = deepcopy(llh_paramset)
    lp = lnprior(theta, paramset=dc_llh_paramset)
    if not np.isfinite(lp):
        return -np.inf
    return lp + triangle_llh(
        theta, args=args, asimov_paramset=dc_asimov_paramset,
        llh_paramset=dc_llh_paramset
    )
