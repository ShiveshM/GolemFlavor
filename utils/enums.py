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
    BASELINE              = 1
    HOLEICE               = 2
    ABSORPTION            = 3
    SCATTERING            = 4
    SCATTERING_ABSORPTION = 5
    STD                   = 6
    STD_HALF1             = 7
    STD_HALF2             = 8
    LONGLIFE              = 9
    DPL                   = 10


class XSCateg(Enum):
    HALF    = 1
    NOM     = 2
    TWICE   = 3
    TWICE01 = 4
    TWICE02 = 5

