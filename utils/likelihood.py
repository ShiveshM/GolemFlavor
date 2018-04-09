# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : April 04, 2018

"""
Likelihood functions for the BSM flavour ratio analysis
"""

from __future__ import absolute_import, division

import argparse
from functools import partial

import numpy as np
from scipy.stats import multivariate_normal

import GolemFitPy as gf

from utils import fr as fr_utils
from utils import gf as gf_utils
from utils.enums import Likelihood, ParamTag
from utils.misc import enum_parse


def gaussian_llh(fr, fr_bf, sigma):
    """Multivariate gaussian likelihood."""
    cov_fr = np.identity(3) * sigma
    return np.log(multivariate_normal.pdf(fr, mean=fr_bf, cov=cov_fr))


def likelihood_argparse(parser):
    parser.add_argument(
        '--likelihood', default='gaussian', type=partial(enum_parse, c=Likelihood),
        choices=Likelihood, help='likelihood contour'
    )


def lnprior(theta, paramset):
    """Priors on theta."""
    ranges = paramset.ranges
    for value, range in zip(theta, ranges):
        if range[0] <= value <= range[1]:
            pass
        else: return -np.inf
    return 0.


def triangle_llh(theta, args, asimov_paramset, mcmc_paramset, fitter):
    """-Log likelihood function for a given theta."""
    if len(theta) != len(mcmc_paramset):
        raise AssertionError(
            'Length of MCMC scan is not the same as the input '
            'params\ntheta={0}\nmcmc_paramset]{1}'.format(theta, mcmc_paramset)
        )
    for idx, param in enumerate(mcmc_paramset):
        param.value = theta[idx]
    hypo_paramset = asimov_paramset
    for param in mcmc_paramset.from_tag(ParamTag.NUISANCE):
        hypo_paramset[param.name].value = param.value

    if args.energy_dependance is EnergyDependance.SPECTRAL:
        bin_centers = np.sqrt(args.binning[:-1]*args.binning[1:])
        bin_width = np.abs(np.diff(args.binning))

    if args.fix_source_ratio:
        if args.energy_dependance is EnergyDependance.MONO:
            source_flux = args.source_ratio
        elif args.energy_dependance is EnergyDependance.SPECTRAL:
	    source_flux = np.array(
		[fr * np.power(bin_centers, SPECTRAL_INDEX)
		 for fr in args.source_ratio]
	    ).T
    else:
        if args.energy_dependance is EnergyDependance.MONO:
            source_flux = fr_utils.angles_to_fr(
                mcmc_paramset.from_tag(ParamTag.SRCANGLES, values=True)
            )
        elif args.energy_dependance is EnergyDependance.SPECTRAL:
            source_flux = np.array(
                [fr * np.power(bin_centers, args.spectral_index)
                 for fr in angles_to_fr(theta[-2:])]
            ).T

    bsm_angles = mcmc_paramset.from_tag(
        [ParamTag.SCALE, ParamTag.MMANGLES], values=True
    )

    if args.energy_dependance is EnergyDependance.MONO:
        u = fr_utils.params_to_BSMu(
            theta      = bsm_angles,
            dim        = args.dimension,
            energy     = args.energy,
            no_bsm     = args.no_bsm,
            fix_mixing = args.fix_mixing,
            fix_scale  = args.fix_scale,
            scale      = args.scale
        )
        fr = fr_utils.u_to_fr(source_flux, u)
    elif args.energy_dependance is EnergyDependance.SPECTRAL:
        mf_perbin = []
        for i_sf, sf_perbin in enumerate(source_flux):
            u = fr_utils.params_to_BSMu(
                theta      = bsm_angles,
                dim        = args.dimension,
                energy     = args.energy,
                no_bsm     = args.no_bsm,
                fix_mixing = args.fix_mixing,
                fix_scale  = args.fix_scale,
                scale      = args.scale
            )
            fr = fr_utils.u_to_fr(sf_perbin, u)
            mf_perbin.append(fr)
        measured_flux = np.array(mf_perbin).T
        intergrated_measured_flux = np.sum(measured_flux * bin_width, axis=1)
        averaged_measured_flux = (1./(args.binning[-1] - args.binning[0])) * \
            intergrated_measured_flux
        fr = averaged_measured_flux / np.sum(averaged_measured_flux)

    for idx, param in enumerate(hypo_paramset.from_tag(ParamTag.BESTFIT)):
        param.value = fr[idx]

    if args.likelihood is Likelihood.FLAT:
        return 1.
    elif args.likelihood is Likelihood.GAUSSIAN:
        fr_bf = args.measured_ratio
        return gaussian_llh(fr, fr_bf, args.sigma_ratio)
    elif args.likelihood is Likelihood.GOLEMFIT:
        return gf_utils.get_llh(fitter, hypo_paramset)
