#! /usr/bin/env python
from __future__ import absolute_import, division

import sys
sys.path.extend(['.', '../'])

import numpy as np

from utils import fr as fr_utils
from utils.enums import MixingScenario

SOURCE = [0, 1, 0]

bsm = True
SCALE = 1E-45
DIMENSION = 6
FIX_MIXING = MixingScenario.T13
ENERGY = 1E6


if len(sys.argv)< 2:
    print sys.argv
    print "Usage: angles_to_fr.py input_filepath."
    exit(1)

infile = sys.argv[1]
outfile = infile[:-4] + '_proc.npy'

d = np.load(infile)

def m_fr(theta):
    if not bsm:
        s_12_2, c_13_4, s_23_2, dcp, m21_2, m3x_2 = theta
        sm_u = fr_utils.angles_to_u((s_12_2, c_13_4, s_23_2, dcp))
        sm_u = np.array(sm_u, dtype=np.complex256)
        return fr_utils.u_to_fr(SOURCE, sm_u)
    elif bsm:
        s_12_2, c_13_4, s_23_2, dcp, m21_2, m3x_2 = theta[:6]
        sm_u = fr_utils.angles_to_u((s_12_2, c_13_4, s_23_2, dcp))
        bsm_u = np.array(
            fr_utils.params_to_BSMu(
                theta[6:], fix_scale=True, scale=SCALE, dim=DIMENSION,
                energy=ENERGY, sm_u=sm_u
            ), dtype=np.complex256
        )
        return fr_utils.u_to_fr(SOURCE, bsm_u)

pd = np.array(map(m_fr, d))
np.save(outfile, pd)
