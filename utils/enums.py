# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : March 17, 2018

"""
Define Enums for the BSM flavour ratio analysis
"""

from enum import Enum


class DataType(Enum):
    REAL    = 1
    FAKE    = 2
    ASMIMOV = 3


class EnergyDependance(Enum):
    MONO     = 1
    SPECTRAL = 2


class Likelihood(Enum):
    FLAT     = 1
    GAUSSIAN = 2
    GOLEMFIT = 3
    GF_FREQ  = 4


class ParamTag(Enum):
    NUISANCE  = 1
    SM_ANGLES = 2
    MMANGLES  = 3
    SCALE     = 4
    SRCANGLES = 5
    BESTFIT   = 6
    NONE      = 7


class PriorsCateg(Enum):
    UNIFORM  = 1
    GAUSSIAN = 2


class MCMCSeedType(Enum):
    UNIFORM  = 1
    GAUSSIAN = 2


class SensitivityCateg(Enum):
    FULL        = 1
    FIXED_ANGLE = 2
    CORR_ANGLE  = 3


class StatCateg(Enum):
    BAYESIAN    = 1
    FREQUENTIST = 2


class SteeringCateg(Enum):
    P2_0 = 1
    P2_1 = 2
