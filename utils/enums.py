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


class Likelihood(Enum):
    GOLEMFIT = 1
    GF_FREQ  = 2


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


class StatCateg(Enum):
    BAYESIAN    = 1
    FREQUENTIST = 2


class SteeringCateg(Enum):
    P2_0 = 1
    P2_1 = 2


class Texture(Enum):
    OEU  = 1
    OET  = 2
    OUT  = 3
    NONE = 4
