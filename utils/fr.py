# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : March 17, 2018

"""
Useful functions for the BSM flavour ratio analysis
"""

from __future__ import absolute_import, division

from functools import partial

import numpy as np

from utils.enums import EnergyDependance, MixingScenario
from utils.misc import enum_parse, parse_bool

import mpmath as mp
mp.mp.dps = 100 # Computation precision

# DTYPE  = np.float128
# CDTYPE = np.complex256
# PI     = np.arccos(DTYPE(-1))
# SQRT   = np.sqrt
# COS    = np.cos
# SIN    = np.sin
# ACOS   = np.arccos
# ASIN   = np.arcsin
# EXP    = np.exp

DTYPE  = mp.mpf
CDTYPE = mp.mpc
PI     = mp.pi
SQRT   = mp.sqrt
COS    = mp.cos
SIN    = mp.sin
ACOS   = mp.acos
ASIN   = mp.asin
EXP    = mp.exp

MASS_EIGENVALUES = [7.40E-23, 2.515E-21]
"""SM mass eigenvalues"""


def determinant(x):
    """Calculate the determininant of a 3x3 matrix.

    Parameters
    ----------
    x : ndarray, shape = (3, 3)

    Returns
    ----------
    float determinant

    Examples
    ----------
    >>> print determinant(
    >>>     [[-1.65238188-0.59549718j,  0.27486548-0.18437467j, -1.35524534-0.38542072j],
    >>>      [-1.07480906+0.29630449j, -0.47808456-0.80316821j, -0.88609356-1.50737308j],
    >>>      [-0.14924144-0.99230446j,  0.49504234+0.63639805j, 2.29258915-0.36537507j]]
    >>> )
    (2.7797571563274688+3.0841795325804848j)

    """
    return (x[0][0] * (x[1][1] * x[2][2] - x[2][1] * x[1][2])
           -x[1][0] * (x[0][1] * x[2][2] - x[2][1] * x[0][2])
           +x[2][0] * (x[0][1] * x[1][2] - x[1][1] * x[0][2]))


def angles_to_fr(src_angles):
    """Convert angular projection of the source flavour ratio back into the
    flavour ratio.

    Parameters
    ----------
    src_angles : list, length = 2
        sin(phi)^4 and cos(psi)^2

    Returns
    ----------
    flavour ratios (nue, numu, nutau)

    Examples
    ----------
    >>> print angles_to_fr((0.3, 0.4))
    (0.38340579025361626, 0.16431676725154978, 0.45227744249483393)

    """
    sphi4, c2psi = map(DTYPE, src_angles)

    psi = (0.5)*ACOS(c2psi)

    sphi2 = SQRT(sphi4)
    cphi2 = 1. - sphi2
    spsi2 = SIN(psi)**2
    cspi2 = 1. - spsi2

    x = sphi2*cspi2
    y = sphi2*spsi2
    z = cphi2
    return x, y, z


def angles_to_u(bsm_angles):
    """Convert angular projection of the mixing matrix elements back into the
    mixing matrix elements.

    Parameters
    ----------
    bsm_angles : list, length = 4
        sin(12)^2, cos(13)^4, sin(23)^2 and deltacp

    Returns
    ----------
    unitary numpy ndarray of shape (3, 3)

    Examples
    ----------
    >>> from fr import angles_to_u
    >>> print angles_to_u((0.2, 0.3, 0.5, 1.5))
    array([[ 0.66195018+0.j        ,  0.33097509+0.j        ,  0.04757188-0.6708311j ],
           [-0.34631487-0.42427084j,  0.61741198-0.21213542j,  0.52331757+0.j        ],
           [ 0.28614067-0.42427084j, -0.64749908-0.21213542j,  0.52331757+0.j        ]])

    """
    s12_2, c13_4, s23_2, dcp = map(DTYPE, bsm_angles)
    dcp = CDTYPE(dcp)

    c12_2 = 1. - s12_2
    c13_2 = SQRT(c13_4)
    s13_2 = 1. - c13_2
    c23_2 = 1. - s23_2

    t12 = ASIN(SQRT(s12_2))
    t13 = ACOS(SQRT(c13_2))
    t23 = ASIN(SQRT(s23_2))

    c12 = COS(t12)
    s12 = SIN(t12)
    c13 = COS(t13)
    s13 = SIN(t13)
    c23 = COS(t23)
    s23 = SIN(t23)

    p1 = np.array([[1   , 0   , 0]                , [0    , c23 , s23] , [0                , -s23 , c23]] , dtype=CDTYPE)
    p2 = np.array([[c13 , 0   , s13*EXP(-1j*dcp)] , [0    , 1   , 0]   , [-s13*EXP(1j*dcp) , 0    , c13]] , dtype=CDTYPE)
    p3 = np.array([[c12 , s12 , 0]                , [-s12 , c12 , 0]   , [0                , 0    , 1]]   , dtype=CDTYPE)

    u = np.dot(np.dot(p1, p2), p3)
    return u


def cardano_eqn(ham):
    """Diagonalise the effective Hamiltonian 3x3 matrix into the form
    h_{eff} = UE_{eff}U^{dagger} using the procedure in PRD91, 052003 (2015).

    Parameters
    ----------
    ham : numpy ndarray of shape (3, 3)
        sin(12)^2, cos(13)^4, sin(23)^2 and deltacp

    Returns
    ----------
    unitary numpy ndarray of shape (3, 3)

    Examples
    ----------
    >>> import numpy as np
    >>> from fr import cardano_eqn
    >>> ham = np.array(
    >>>     [[ 0.66195018+0.j        ,  0.33097509+0.j        ,  0.04757188-0.6708311j ],
    >>>      [-0.34631487-0.42427084j,  0.61741198-0.21213542j,  0.52331757+0.j        ],
    >>>      [ 0.28614067-0.42427084j, -0.64749908-0.21213542j,  0.52331757+0.j        ]]
    >>> )
    >>> print cardano_eqn(ham)
    array([[-0.11143379-0.58863683j, -0.09067747-0.48219068j, 0.34276625-0.08686465j],
           [ 0.14835519+0.47511473j, -0.18299305+0.40777481j, 0.31906300+0.82514223j],
           [-0.62298966+0.07231745j, -0.61407815-0.42709603j, 0.03660313+0.30160428j]])

    """
    if np.shape(ham) != (3, 3):
        raise ValueError(
            'Input matrix should be a square and dimension 3, '
            'got\n{0}'.format(ham)
        )

    a = -np.trace(ham)
    b = DTYPE(1)/2 * ((np.trace(ham))**DTYPE(2) - np.trace(np.dot(ham, ham)))
    c = -determinant(ham)

    Q = (DTYPE(1)/9) * (a**DTYPE(2) - DTYPE(3)*b)
    R = (DTYPE(1)/54) * (DTYPE(2)*a**DTYPE(3) - DTYPE(9)*a*b + DTYPE(27)*c)
    theta = ACOS(R / SQRT(Q**DTYPE(3)))

    E1 = -DTYPE(2) * SQRT(Q) * COS(theta/DTYPE(3)) - (DTYPE(1)/3)*a
    E2 = -DTYPE(2) * SQRT(Q) * COS((theta - DTYPE(2)*PI)/DTYPE(3)) - (DTYPE(1)/3)*a
    E3 = -DTYPE(2) * SQRT(Q) * COS((theta + DTYPE(2)*PI)/DTYPE(3)) - (DTYPE(1)/3)*a

    A1 = ham[1][2] * (ham[0][0] - E1) - ham[1][0]*ham[0][2]
    A2 = ham[1][2] * (ham[0][0] - E2) - ham[1][0]*ham[0][2]
    A3 = ham[1][2] * (ham[0][0] - E3) - ham[1][0]*ham[0][2]

    B1 = ham[2][0] * (ham[1][1] - E1) - ham[2][1]*ham[1][0]
    B2 = ham[2][0] * (ham[1][1] - E2) - ham[2][1]*ham[1][0]
    B3 = ham[2][0] * (ham[1][1] - E3) - ham[2][1]*ham[1][0]

    C1 = ham[1][0] * (ham[2][2] - E1) - ham[1][2]*ham[2][0]
    C2 = ham[1][0] * (ham[2][2] - E2) - ham[1][2]*ham[2][0]
    C3 = ham[1][0] * (ham[2][2] - E3) - ham[1][2]*ham[2][0]

    N1 = SQRT(np.abs(A1*B1)**2 + np.abs(A1*C1)**2 + np.abs(B1*C1)**2)
    N2 = SQRT(np.abs(A2*B2)**2 + np.abs(A2*C2)**2 + np.abs(B2*C2)**2)
    N3 = SQRT(np.abs(A3*B3)**2 + np.abs(A3*C3)**2 + np.abs(B3*C3)**2)

    mm = np.array([
        [np.conjugate(B1)*C1 / N1, np.conjugate(B2)*C2 / N2, np.conjugate(B3)*C3 / N3],
        [A1*C1 / N1, A2*C2 / N2, A3*C3 / N3],
        [A1*B1 / N1, A2*B2 / N2, A3*B3 / N3]
    ])
    return mm


def normalise_fr(fr):
    """Normalise an input flavour combination to a flavour ratio.

    Parameters
    ----------
    fr : list, length = 3
        flavour combination

    Returns
    ----------
    numpy ndarray flavour ratio

    Examples
    ----------
    >>> from fr import normalise_fr
    >>> print normalise_fr((1, 2, 3))
    array([ 0.16666667,  0.33333333,  0.5       ])

    """
    return np.array(fr) / float(np.sum(fr))


def estimate_scale(args):
    """Estimate the scale at which new physics will enter."""
    try: m_eign = args.m3x_2
    except: m_eign = MASS_EIGENVALUES[1]
    if hasattr(args, 'scale'):
        if args.scale != 0:
            scale = args.scale
            scale_region = (scale/args.scale_region, scale*args.scale_region)
            return scale, scale_region
    if args.energy_dependance is EnergyDependance.MONO:
        scale = np.power(
            10, np.round(np.log10(m_eign/args.energy)) - \
            np.log10(args.energy**(args.dimension-3))
        )
        scale_region = (scale/args.scale_region, scale*args.scale_region)
    elif args.energy_dependance is EnergyDependance.SPECTRAL:
        lower_s = (m_eign/args.binning[-1]) / (args.binning[-1]**(args.dimension-3))
        upper_s = (m_eign/args.binning[0]) / (args.binning[0]**(args.dimension-3))
        scale = np.power(10, np.average(np.log10([lower_s, upper_s])))
        diff = upper_s / lower_s
        scale_region = (lower_s/np.power(10, args.dimension), upper_s*diff*np.power(10, args.dimension))
        scale_region = [np.power(10, np.round(np.log10(x))) for x in scale_region]
    return scale, scale_region


def fr_argparse(parser):
    parser.add_argument(
        '--measured-ratio', type=float, nargs=3, default=[1, 1, 1],
        help='Set the central value for the measured flavour ratio at IceCube'
    )
    parser.add_argument(
        '--source-ratio', type=float, nargs=3, default=[2, 1, 0],
        help='Set the source flavour ratio for the case when you want to fix it'
    )
    parser.add_argument(
        '--no-bsm', type=parse_bool, default='False',
        help='Turn off BSM terms'
    )
    parser.add_argument(
        '--dimension', type=int, default=3,
        help='Set the new physics dimension to consider'
    )
    parser.add_argument(
        '--energy-dependance', default='spectral', type=partial(enum_parse, c=EnergyDependance),
        choices=EnergyDependance,
        help='Type of energy dependance to use in the BSM fit'
    )
    parser.add_argument(
        '--spectral-index', default=-2, type=float,
        help='Spectral index for spectral energy dependance'
    )
    parser.add_argument(
        '--fold-index', default='True', type=parse_bool,
        help='Fold in the spectral index when using GolemFit'
    )
    parser.add_argument(
        '--binning', default=[6e4, 1e7, 20], type=float, nargs=3,
        help='Binning for spectral energy dependance'
    )
    parser.add_argument(
        '--fix-source-ratio', type=parse_bool, default='False',
        help='Fix the source flavour ratio'
    )
    parser.add_argument(
        '--fix-mixing', type=partial(enum_parse, c=MixingScenario),
        default='None', choices=MixingScenario,
        help='Fix all mixing parameters to choice of maximal mixing'
    )
    parser.add_argument(
        '--fix-mixing-almost', type=parse_bool, default='False',
        help='Fix all mixing parameters except s23'
    )
    parser.add_argument(
        '--fix-scale', type=parse_bool, default='False',
        help='Fix the new physics scale'
    )
    parser.add_argument(
        '--scale', type=float, default=0,
        help='Set the new physics scale'
    )
    parser.add_argument(
        '--scale-region', type=float, default=1e10,
        help='Set the size of the box to scan for the scale'
    )
    parser.add_argument(
        '--energy', type=float, default=1000,
        help='Set the energy scale'
    )


def fr_to_angles(ratios):
    """Convert from flavour ratio into the angular projection of the flavour
    ratios.

    Parameters
    ----------
    TODO(shivesh)
    """
    fr0, fr1, fr2 = normalise_fr(ratios)

    cphi2 = fr2
    sphi2 = (1.0 - cphi2)

    if sphi2 == 0.:
        return (0., 0.)
    else:
        cpsi2 = fr0 / sphi2

    sphi4 = sphi2**2
    c2psi = COS(ACOS(SQRT(cpsi2))*2)

    return map(float, (sphi4, c2psi))


NUFIT_U = angles_to_u((0.307, (1-0.02195)**2, 0.565, 3.97935))
"""NuFIT mixing matrix (s_12^2, c_13^4, s_23^2, dcp)"""


def params_to_BSMu(theta, dim, energy, mass_eigenvalues=MASS_EIGENVALUES,
                   sm_u=NUFIT_U, no_bsm=False, fix_mixing=MixingScenario.NONE,
                   fix_mixing_almost=False, fix_scale=False, scale=None,
                   check_uni=True, epsilon=1e-7):
    """Construct the BSM mixing matrix from the BSM parameters.

    Parameters
    ----------
    theta : list, length > 3
        BSM parameters

    dim : int
        Dimension of BSM physics

    energy : float
        Energy in GeV

    mass_eigenvalues : list, length = 2
        SM mass eigenvalues

    sm_u : numpy ndarray, dimension 3
        SM mixing matrix

    no_bsm : bool
        Turn off BSM behaviour

    fix_mixing : MixingScenario
        Fix the BSM mixing angles

    fix_mixing_almost : bool
        Fix the BSM mixing angles except one

    fix_scale : bool
        Fix the BSM scale

    scale : float
        Used with fix_scale - scale at which to fix

    check_uni : bool
        Check the resulting BSM mixing matrix is unitary

    Returns
    ----------
    unitary numpy ndarray of shape (3, 3)

    Examples
    ----------
    >>> from fr import params_to_BSMu
    >>> print params_to_BSMu((0.2, 0.3, 0.5, 1.5, -20), dim=3, energy=1000)
    array([[ 0.18658169 -6.34190523e-01j, -0.26460391 +2.01884200e-01j, 0.67247096 -9.86808417e-07j],
           [-0.50419832 +2.14420570e-01j, -0.36013768 +5.44254868e-01j, 0.03700961 +5.22039894e-01j],
           [-0.32561308 -3.95946524e-01j,  0.64294909 -2.23453580e-01j, 0.03700830 +5.22032403e-01j]])

    """
    if np.shape(sm_u) != (3, 3):
        raise ValueError(
            'Input matrix should be a square and dimension 3, '
            'got\n{0}'.format(sm_u)
        )

    if fix_mixing is not MixingScenario.NONE and fix_mixing_almost:
        raise NotImplementedError(
            '--fix-mixing and --fix-mixing-almost cannot be used together'
        )

    if not isinstance(theta, (list, tuple)):
        theta = [theta]

    if fix_mixing is MixingScenario.T12:
        s12_2, c13_4, s23_2, dcp, sc2 = 0.5, 1.0, 0.0, 0., theta
    elif fix_mixing is MixingScenario.T13:
        s12_2, c13_4, s23_2, dcp, sc2 = 0.0, 0.25, 0.0, 0., theta
    elif fix_mixing is MixingScenario.T23:
        s12_2, c13_4, s23_2, dcp, sc2 = 0.0, 1.0, 0.5, 0., theta
    elif fix_mixing_almost:
        s12_2, c13_4, dcp = 0.5, 1.0-1E-6, 0.
        s23_2, sc2 = theta
    elif fix_scale:
        s12_2, c13_4, s23_2, dcp = theta
        sc2 = scale
    else:
        s12_2, c13_4, s23_2, dcp, sc2 = theta

    if len(theta) != 0:
        sc2 = np.power(10., sc2)
    else:
        sc2 = scale
    sc1 = sc2 / 100.

    mass_matrix = np.array(
        [[0, 0, 0], [0, mass_eigenvalues[0], 0], [0, 0, mass_eigenvalues[1]]]
    )
    sm_ham = (1./(2*energy))*np.dot(sm_u, np.dot(mass_matrix, sm_u.conj().T))
    if no_bsm:
        eg_vector = cardano_eqn(sm_ham)
    else:
        new_physics_u = angles_to_u((s12_2, c13_4, s23_2, dcp))
        scale_matrix = np.array(
            [[0, 0, 0], [0, sc1, 0], [0, 0, sc2]]
        )
        bsm_term = (energy**(dim-3)) * np.dot(new_physics_u, np.dot(scale_matrix, new_physics_u.conj().T))
        bsm_ham = sm_ham + bsm_term
        eg_vector = cardano_eqn(bsm_ham)

    if check_uni:
        test_unitarity(eg_vector, rse=True, epsilon=epsilon)
    return eg_vector


def test_unitarity(x, prnt=False, rse=False, epsilon=None):
    """Test the unitarity of a matrix.

    Parameters
    ----------
    x : numpy ndarray
        Matrix to evaluate

    prnt : bool
        Print the result

    rse : bool
        Raise Assertion if matrix is not unitary

    Returns
    ----------
    numpy ndarray

    Examples
    ----------
    >>> from fr import test_unitarity
    >>> x = np.identity(3)
    >>> print test_unitarity(x)
    array([[ 1.,  0.,  0.],
           [ 0.,  1.,  0.],
           [ 0.,  0.,  1.]])

    """
    f = np.abs(np.dot(x, x.conj().T), dtype=DTYPE)
    if prnt:
        print 'Unitarity test:\n{0}'.format(f)
    if rse:
        if not np.abs(np.trace(f) - 3.) < epsilon or \
           not np.abs(np.sum(f) - 3.) < epsilon:
            raise AssertionError(
                'Matrix is not unitary!\nx\n{0}\ntest '
                'u\n{1}'.format(x, f)
            )
    return f


def u_to_fr(source_fr, matrix):
    """Compute the observed flavour ratio assuming decoherence.

    Parameters
    ----------
    source_fr : list, length = 3
        Source flavour ratio components

    matrix : numpy ndarray, dimension 3
        Mixing matrix

    Returns
    ----------
    Measured flavour ratio

    Examples
    ----------
    >>> from fr import params_to_BSMu, u_to_fr
    >>> print u_to_fr((1, 2, 0), params_to_BSMu((0.2, 0.3, 0.5, 1.5, -20), 3, 1000))
        array([ 0.33740075,  0.33176584,  0.33083341])

    """
    try:
        composition = np.einsum(
            'ai, bi, a -> b', np.abs(matrix)**2, np.abs(matrix)**2, source_fr,
        )
    except:
        matrix = np.array(matrix, dtype=np.complex256)
        composition = np.einsum(
            'ai, bi, a -> b', np.abs(matrix)**2, np.abs(matrix)**2, source_fr,
        )
        pass

    ratio = composition / np.sum(source_fr)
    return ratio
