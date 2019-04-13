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

from utils import fr as fr_utils
from utils import gf as gf_utils
from utils.enums import Likelihood, ParamTag, PriorsCateg, StatCateg
from utils.misc import enum_parse, gen_identifier, parse_bool


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
    for idx, param in enumerate(llh_paramset):
        param.value = theta[idx]
    hypo_paramset = asimov_paramset
    for param in llh_paramset.from_tag(ParamTag.NUISANCE):
        hypo_paramset[param.name].value = param.value

    bin_centers = np.sqrt(args.binning[:-1]*args.binning[1:])
    bin_width = np.abs(np.diff(args.binning))
    spectral_index = -hypo_paramset['astroDeltaGamma'].value

    source_flux = np.array(
        [fr * np.power(bin_centers, spectral_index)
         for fr in args.source_ratio]
    ).T

    bsm_angles = llh_paramset.from_tag(
        [ParamTag.SCALE, ParamTag.MMANGLES], values=True
    )

    m_eig_names = ['m21_2', 'm3x_2']
    ma_names = ['s_12_2', 'c_13_4', 's_23_2', 'dcp']

    if set(m_eig_names+ma_names).issubset(set(llh_paramset.names)):
        mass_eigenvalues = [x.value for x in llh_paramset if x.name in m_eig_names]
        sm_u = fr_utils.angles_to_u(
            [x.value for x in llh_paramset if x.name in ma_names]
        )
    else:
        mass_eigenvalues = fr_utils.MASS_EIGENVALUES
        sm_u = fr_utils.NUFIT_U

    if args.no_bsm:
        fr = fr_utils.u_to_fr(source_flux, np.array(sm_u, dtype=np.complex256))
    else:
        mf_perbin = []
        for i_sf, sf_perbin in enumerate(source_flux):
            u = fr_utils.params_to_BSMu(
                theta             = bsm_angles,
                dim               = args.dimension,
                energy            = bin_centers[i_sf],
                mass_eigenvalues  = mass_eigenvalues,
                sm_u              = sm_u,
                no_bsm            = args.no_bsm,
                texture           = args.texture,
            )
            fr = fr_utils.u_to_fr(sf_perbin, u)
            mf_perbin.append(fr)
        measured_flux = np.array(mf_perbin).T
        intergrated_measured_flux = np.sum(measured_flux * bin_width, axis=1)
        averaged_measured_flux = (1./(args.binning[-1] - args.binning[0])) * \
            intergrated_measured_flux
        fr = averaged_measured_flux / np.sum(averaged_measured_flux)

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
