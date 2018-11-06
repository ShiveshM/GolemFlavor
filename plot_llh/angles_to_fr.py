#! /usr/bin/env python
from __future__ import absolute_import, division

import sys
sys.path.extend(['.', '../'])

import numpy as np

from utils import fr as fr_utils

SOURCE = [1, 2, 0]

if len(sys.argv)< 2:
    print sys.argv
    print "Usage: angles_to_fr.py input_filepath."
    exit(1)

infile = sys.argv[1]
outfile = infile[:-4] + '_proc.npy'

d = np.load(infile)

def m_fr(theta):
    s_12_2, c_13_4, s_23_2, dcp, m21_2, m3x_2 = theta
    sm_u = np.array(
        fr_utils.angles_to_u((s_12_2, c_13_4, s_23_2, dcp)), dtype=np.complex256
    )
    return fr_utils.u_to_fr(SOURCE, sm_u)

pd = np.array(map(m_fr, d))
np.save(outfile, pd)
