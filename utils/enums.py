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


class FitCateg(Enum):
    HESESPL  = 1
    HESEDPL  = 2
    ZPSPL    = 3
    ZPDPL    = 4
    NUNUBAR2 = 5


class FitFlagCateg(Enum):
    DEFAULT = 1
    XS      = 2


class FitPriorsCateg(Enum):
    DEFAULT = 1
    XS      = 2


class Likelihood(Enum):
    FLAT     = 1
    GAUSSIAN = 2
    GOLEMFIT = 3


class ParamTag(Enum):
    NUISANCE  = 1
    MMANGLES  = 2
    SCALE     = 3
    SRCANGLES = 4
    BESTFIT   = 5
    NONE      = 6


class Priors(Enum):
    UNIFORM    = 1
    LOGUNIFORM = 2
    JEFFREYS   = 3


class SteeringCateg(Enum):
    P2_0                  = 1
    P2_1                  = 2
    # TODO(shivesh): fix this "P2_-1"
    P2__1                 = 3
    P2__3                 = 4
    P2_0_HALF1            = 5
    P2_0_HALF2            = 6
    ABSORPTION            = 7
    SCATTERING            = 8
    SCATTERING_ABSORPTION = 9
    LONGLIFE              = 10
    DPL                   = 11

class XSCateg(Enum):
    HALF    = 1
    NOM     = 2
    TWICE   = 3
    TWICE01 = 4
    TWICE02 = 5

