#! /usr/bin/env python
from __future__ import absolute_import, division

import sys
sys.path.extend(['.', '..'])

import numpy as np
import tqdm

from utils import fr as fr_utils
from utils import misc as misc_utils
from utils.enums import MixingScenario

binning = np.logspace(np.log10(6e4), np.log10(1e7), 21)
dimension = 6
source = [0, 1, 0]
scenario = MixingScenario.T13

def get_fr(theta, source, binning, dimension, scenario):
    sm_mixings = theta[:6]
    nuisance   = theta[6:11]
    scale      = np.power(10., theta[-1])

    index = -nuisance[-1]

    bin_centers = np.sqrt(binning[:-1]*binning[1:])
    bin_width = np.abs(np.diff(binning))

    # TODO(shivesh): test with astroNorm
    source_flux = np.array(
        [fr * np.power(bin_centers, index) for fr in source]
    ).T

    mass_eigenvalues = sm_mixings[-2:]
    sm_u = fr_utils.angles_to_u(sm_mixings[:-2])

    mf_perbin = []
    for i_sf, sf_perbin in enumerate(source_flux):
        u = fr_utils.params_to_BSMu(
            theta             = [],
            dim               = dimension,
            energy            = bin_centers[i_sf],
            mass_eigenvalues  = mass_eigenvalues,
            sm_u              = sm_u,
            no_bsm            = False,
            fix_mixing        = scenario,
            fix_mixing_almost = False,
            fix_scale         = True,
            scale             = scale
        )
        fr = fr_utils.u_to_fr(sf_perbin, u)
        mf_perbin.append(fr)
    measured_flux = np.array(mf_perbin).T
    intergrated_measured_flux = np.sum(measured_flux * bin_width, axis=1)
    averaged_measured_flux = (1./(binning[-1] - binning[0])) * \
        intergrated_measured_flux
    fr = averaged_measured_flux / np.sum(averaged_measured_flux)

    return map(float, fr)

if len(sys.argv)< 2:
    print sys.argv
    print "Usage: calc_fr_MCMC.py input_filepath."
    exit(1)

infile = sys.argv[1]
outfile = infile[:-4] + '_proc.npy'

d = np.load(infile)

p = []
for x in tqdm.tqdm(d, total=len(d)):
    p.append(get_fr(x, source, binning, dimension, scenario))
p = np.array(p)

np.save(outfile, p)
