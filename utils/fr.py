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

from utils.enums import EnergyDependance
from utils.misc import DTYPE, CDTYPE, PI
from utils.misc import enum_parse, parse_bool


MASS_EIGENVALUES = [7.40E-23, 2.515E-21]
"""SM mass eigenvalues"""


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
    sphi4, c2psi = src_angles

    psi = (0.5)*np.arccos(c2psi, dtype=DTYPE)

    sphi2 = np.sqrt(sphi4, dtype=DTYPE)
    cphi2 = 1. - sphi2
    spsi2 = np.sin(psi, dtype=DTYPE)**2
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
    s12_2, c13_4, s23_2, dcp = bsm_angles
    dcp = CDTYPE(dcp)

    c12_2 = 1. - s12_2
    c13_2 = np.sqrt(c13_4, dtype=DTYPE)
    s13_2 = 1. - c13_2
    c23_2 = 1. - s23_2

    t12 = np.arcsin(np.sqrt(s12_2, dtype=DTYPE))
    t13 = np.arccos(np.sqrt(c13_2, dtype=DTYPE))
    t23 = np.arcsin(np.sqrt(s23_2, dtype=DTYPE))

    c12 = np.cos(t12)
    s12 = np.sin(t12)
    c13 = np.cos(t13)
    s13 = np.sin(t13)
    c23 = np.cos(t23)
    s23 = np.sin(t23)

    p1 = np.array([[1   , 0   , 0]                   , [0    , c23 , s23] , [0                   , -s23 , c23]] , dtype=CDTYPE)
    p2 = np.array([[c13 , 0   , s13*np.exp(-1j*dcp)] , [0    , 1   , 0]   , [-s13*np.exp(1j*dcp) , 0    , c13]] , dtype=CDTYPE)
    p3 = np.array([[c12 , s12 , 0]                   , [-s12 , c12 , 0]   , [0                   , 0    , 1]]   , dtype=CDTYPE)

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
    b = (0.5) * ((np.trace(ham))**2 - np.trace(np.dot(ham, ham)))
    c = -np.linalg.det(ham)

    Q = (1/9.) * (a**2 - 3*b)
    R = (1/54.) * (2*a**3 - 9*a*b + 27*c)
    theta = np.arccos(R / np.sqrt(Q**3))

    E1 = -2 * np.sqrt(Q) * np.cos(theta/3.) - (1/3.)*a
    E2 = -2 * np.sqrt(Q) * np.cos((theta - 2.*PI)/3.) - (1/3.)*a
    E3 = -2 * np.sqrt(Q) * np.cos((theta + 2.*PI)/3.) - (1/3.)*a

    A1 = ham[1][2] * (ham[0][0] - E1) - ham[1][0]*ham[0][2]
    A2 = ham[1][2] * (ham[0][0] - E2) - ham[1][0]*ham[0][2]
    A3 = ham[1][2] * (ham[0][0] - E3) - ham[1][0]*ham[0][2]

    B1 = ham[2][0] * (ham[1][1] - E1) - ham[2][1]*ham[1][0]
    B2 = ham[2][0] * (ham[1][1] - E2) - ham[2][1]*ham[1][0]
    B3 = ham[2][0] * (ham[1][1] - E3) - ham[2][1]*ham[1][0]

    C1 = ham[1][0] * (ham[2][2] - E1) - ham[1][2]*ham[2][0]
    C2 = ham[1][0] * (ham[2][2] - E2) - ham[1][2]*ham[2][0]
    C3 = ham[1][0] * (ham[2][2] - E3) - ham[1][2]*ham[2][0]

    N1 = np.sqrt(np.abs(A1*B1, dtype=DTYPE)**2 + np.abs(A1*C1, dtype=DTYPE)**2 + np.abs(B1*C1, dtype=DTYPE)**2)
    N2 = np.sqrt(np.abs(A2*B2, dtype=DTYPE)**2 + np.abs(A2*C2, dtype=DTYPE)**2 + np.abs(B2*C2, dtype=DTYPE)**2)
    N3 = np.sqrt(np.abs(A3*B3, dtype=DTYPE)**2 + np.abs(A3*C3, dtype=DTYPE)**2 + np.abs(B3*C3, dtype=DTYPE)**2)

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
        scale_region = (lower_s/diff, upper_s*diff)
        scale_region = [np.power(10, np.round(np.log10(x))) for x in scale_region]
    return scale, scale_region


def fr_argparse(parser):
    parser.add_argument(
        '--measured-ratio', type=int, nargs=3, default=[1, 1, 1],
        help='Set the central value for the measured flavour ratio at IceCube'
    )
    parser.add_argument(
        '--source-ratio', type=int, nargs=3, default=[2, 1, 0],
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
        '--spectral-index', default=-2, type=int,
        help='Spectral index for spectral energy dependance'
    )
    parser.add_argument(
        '--fold-index', default='True', type=parse_bool,
        help='Fold in the spectral index when using GolemFit'
    )
    parser.add_argument(
        '--binning', default=[1e4, 1e7, 5], type=float, nargs=3,
        help='Binning for spectral energy dependance'
    )
    parser.add_argument(
        '--fix-source-ratio', type=parse_bool, default='False',
        help='Fix the source flavour ratio'
    )
    parser.add_argument(
        '--fix-mixing', type=parse_bool, default='False',
        help='Fix all mixing parameters to bi-maximal mixing'
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
        '--scale', type=float, default=1e-30,
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
    sphi4 = (fr2 - 1.0)**2
    if (fr2 - 1.0) == 0:
        c2psi = 0
    else:
        c2psi = (fr1*2.0 + fr2 - 1.0) * (fr2 - 1.0)
    return sphi4, c2psi


NUFIT_U = angles_to_u((0.307, (1-0.02195)**2, 0.565, 3.97935))
"""NuFIT mixing matrix (s_12^2, c_13^4, s_23^2, dcp)"""


def params_to_BSMu(theta, dim, energy, mass_eigenvalues=MASS_EIGENVALUES,
                   sm_u=NUFIT_U, no_bsm=False, fix_mixing=False,
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

    fix_mixing : bool
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

    if fix_mixing and fix_mixing_almost:
        raise NotImplementedError(
            '--fix-mixing and --fix-mixing-almost cannot be used together'
        )

    if fix_mixing:
        s12_2, c13_4, s23_2, dcp, sc2 = 0.5, 1.0-1E-6, 0.5, 0., theta
    elif fix_mixing_almost:
        s12_2, c13_4, dcp = 0.5, 1.0-1E-6, 0.
        s23_2, sc2 = theta
    elif fix_scale:
        s12_2, c13_4, s23_2, dcp = theta
        sc2 = np.log10(scale)
    else:
        s12_2, c13_4, s23_2, dcp, sc2 = theta
    sc2 = np.power(10., sc2)
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
        tu = test_unitarity(eg_vector)
        if not np.abs(np.trace(tu) - 3.) < epsilon or \
           not np.abs(np.sum(tu) - 3.) < epsilon:
            raise AssertionError(
                'Matrix is not unitary!\neg_vector\n{0}\ntest '
                'u\n{1}'.format(eg_vector, tu)
            )
    return eg_vector


def test_unitarity(x, prnt=False):
    """Test the unitarity of a matrix.

    Parameters
    ----------
    x : numpy ndarray
        Matrix to evaluate

    prnt : bool
        Print the result

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
    composition = np.einsum(
        'ai, bi, a -> b', np.abs(matrix)**2, np.abs(matrix)**2, source_fr,
    )
    ratio = composition / np.sum(source_fr)
    return ratio
