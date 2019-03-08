# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : March 17, 2018

"""
Define Enums for the BSM flavour ratio analysis
"""

from enum import Enum


def str_enum(x):
    return '{0}'.format(str(x).split('.')[-1])


class DataType(Enum):
    REAL        = 1
    ASIMOV      = 2
    REALISATION = 3


class EnergyDependance(Enum):
    MONO     = 1
    SPECTRAL = 2


class Likelihood(Enum):
    FLAT     = 1
    GAUSSIAN = 2
    GOLEMFIT = 3
    GF_FREQ  = 4


class MixingScenario(Enum):
    T12  = 1
    T13  = 2
    T23  = 3
    NONE = 4


class ParamTag(Enum):
    NUISANCE  = 1
    SM_ANGLES = 2
    MMANGLES  = 3
    SCALE     = 4
    SRCANGLES = 5
    BESTFIT   = 6
    NONE      = 7


class PriorsCateg(Enum):
    UNIFORM      = 1
    GAUSSIAN     = 2
    LIMITEDGAUSS = 3


class MCMCSeedType(Enum):
    UNIFORM  = 1
    GAUSSIAN = 2


class SensitivityCateg(Enum):
    FULL            = 1
    FIXED_ANGLE     = 2
    CORR_ANGLE      = 3
    FIXED_ONE_ANGLE = 4
    CORR_ONE_ANGLE  = 5


class StatCateg(Enum):
    BAYESIAN    = 1
    FREQUENTIST = 2


class SteeringCateg(Enum):
    P2_0 = 1
    P2_1 = 2
