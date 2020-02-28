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

from golemflavor.enums import ParamTag, Texture
from golemflavor.misc import enum_parse, parse_bool

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
"""SM mass eigenvalues."""

SCALE_BOUNDARIES = {
    3 : (-32, -20),
    4 : (-40, -24),
    5 : (-48, -27),
    6 : (-56, -30),
    7 : (-64, -33),
    8 : (-72, -36)
}
"""Boundaries to scan the NP scale for each dimension."""


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

    x = float(abs(sphi2*cspi2))
    y = float(abs(sphi2*spsi2))
    z = float(abs(cphi2))
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


def fr_argparse(parser):
    parser.add_argument(
        '--injected-ratio', type=float, nargs=3, required=False,
        help='Injected ratio if not using data'
    )
    parser.add_argument(
        '--source-ratio', type=float, nargs=3, default=[1, 2, 0],
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
        '--texture', type=partial(enum_parse, c=Texture),
        default='none', choices=Texture, help='Set the BSM mixing texture'
    )
    parser.add_argument(
        '--binning', default=[6e4, 1e7, 20], type=float, nargs=3,
        help='Binning for spectral energy dependance'
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


def params_to_BSMu(bsm_angles, dim, energy, mass_eigenvalues=MASS_EIGENVALUES,
                   sm_u=NUFIT_U, no_bsm=False, texture=Texture.NONE,
                   check_uni=True, epsilon=1e-7):
    """Construct the BSM mixing matrix from the BSM parameters.

    Parameters
    ----------
    bsm_angles : list, length > 3
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

    texture : Texture
        BSM mixing texture

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

    if not isinstance(bsm_angles, (list, tuple)):
        bsm_angles = [bsm_angles]

    z = 0.+1e-9
    if texture is Texture.OEU:
        np_s12_2, np_c13_4, np_s23_2, np_dcp, sc2 = 0.5, 1.0, z, z, bsm_angles
    elif texture is Texture.OET:
        np_s12_2, np_c13_4, np_s23_2, np_dcp, sc2 = z, 0.25,  z, z, bsm_angles
    elif texture is Texture.OUT:
        np_s12_2, np_c13_4, np_s23_2, np_dcp, sc2 = z, 1.0, 0.5, z, bsm_angles
    else:
        np_s12_2, np_c13_4, np_s23_2, np_dcp, sc2 = bsm_angles

    sc2 = np.power(10., sc2)
    sc1 = sc2 / 100.

    mass_matrix = np.array(
        [[0, 0, 0], [0, mass_eigenvalues[0], 0], [0, 0, mass_eigenvalues[1]]]
    )
    sm_ham = (1./(2*energy))*np.dot(sm_u, np.dot(mass_matrix, sm_u.conj().T))
    if no_bsm:
        eg_vector = cardano_eqn(sm_ham)
    else:
        NP_U = angles_to_u((np_s12_2, np_c13_4, np_s23_2, np_dcp))
        SC_U = np.array(
            [[0, 0, 0], [0, sc1, 0], [0, 0, sc2]]
        )
        bsm_term = (energy**(dim-3)) * np.dot(NP_U, np.dot(SC_U, NP_U.conj().T))
        bsm_ham = sm_ham + bsm_term
        eg_vector = cardano_eqn(bsm_ham)

    if check_uni:
        test_unitarity(eg_vector, rse=True, epsilon=epsilon)
    return eg_vector


def flux_averaged_BSMu(theta, args, spectral_index, llh_paramset):
    if len(theta) != len(llh_paramset):
        raise AssertionError(
            'Length of MCMC scan is not the same as the input '
            'params\ntheta={0}\nparamset]{1}'.format(theta, llh_paramset)
        )

    for idx, param in enumerate(llh_paramset):
        param.value = theta[idx]

    bin_centers = np.sqrt(args.binning[:-1]*args.binning[1:])
    bin_width = np.abs(np.diff(args.binning))

    source_flux = np.array(
        [fr * np.power(bin_centers, spectral_index)
         for fr in args.source_ratio]
    ).T

    bsm_angles = llh_paramset.from_tag(
        [ParamTag.SCALE, ParamTag.MMANGLES], values=True
    )

    m_eig_names = ['m21_2', 'm3x_2']
    ma_names = ['s_12_2', 'c_13_4', 's_23_2', 'dcp']

    if set(m_eig_names+ma_names).issubset(set(llh_paramset.names)):
        mass_eigenvalues = [x.value for x in llh_paramset if x.name in m_eig_names]
        sm_u = angles_to_u(
            [x.value for x in llh_paramset if x.name in ma_names]
        )
    else:
        mass_eigenvalues = MASS_EIGENVALUES
        sm_u = NUFIT_U

    if args.no_bsm:
        fr = u_to_fr(source_flux, np.array(sm_u, dtype=np.complex256))
    else:
        mf_perbin = []
        for i_sf, sf_perbin in enumerate(source_flux):
            u = params_to_BSMu(
                bsm_angles        = bsm_angles,
                dim               = args.dimension,
                energy            = bin_centers[i_sf],
                mass_eigenvalues  = mass_eigenvalues,
                sm_u              = sm_u,
                no_bsm            = args.no_bsm,
                texture           = args.texture,
            )
            fr = u_to_fr(sf_perbin, u)
            mf_perbin.append(fr)
        measured_flux = np.array(mf_perbin).T
        intergrated_measured_flux = np.sum(measured_flux * bin_width, axis=1)
        averaged_measured_flux = (1./(args.binning[-1] - args.binning[0])) * \
            intergrated_measured_flux
        fr = averaged_measured_flux / np.sum(averaged_measured_flux)
    return fr


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
