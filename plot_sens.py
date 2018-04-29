#! /usr/bin/env python
# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : April 28, 2018

"""
HESE BSM flavour ratio analysis plotting script
"""

from __future__ import absolute_import, division

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import rc
from matplotlib import pyplot as plt
from matplotlib.offsetbox import AnchoredText

rc('text', usetex=False)
rc('font', **{'family':'serif', 'serif':['Computer Modern'], 'size':18})

FRS = [
    (1, 1, 1, 1, 2, 0),
    (1, 1, 1, 0, 1, 0),
    # (1, 1, 1, 1, 0, 0),
    # (1, 1, 1, 0, 0, 1),
]

DIMENSION = [3, 6]
