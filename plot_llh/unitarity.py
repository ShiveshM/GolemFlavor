# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : Sept 11, 2019

"""
Utility methods for calculation of unitarity bounds.
Calculation follows from DOI 10.1103/PhysRevD.98.123023
"Unitary bounds of astrophysical neutrinos" by M. Ahlers, M. Bustamante, S. Mu
"""

from __future__ import absolute_import, division

from copy import deepcopy

import numpy as np

def Sp(x, y, z):
    """S'."""
    if x == 0: return -np.inf
    if np.power(x, 2) < (np.power(y - z, 2) / 9.):
        return -np.inf
    num = np.power(3*x + y + z, 2) - (4 * y * z)
    den = 24 * x
    return num / den

S1 = lambda x, y, z: (x + y + z) / 3.
S2 = lambda x, y, z: [x/2., y/2., z/2.]
S3 = lambda x, y, z: [Sp(x, y, z), Sp(y, z, x), Sp(z, x, y)]

B_v = np.vectorize(
    lambda x, y, z: np.max([0] + [S1(x, y, z)] + S2(x, y, z) + S3(x, y, z))
)

def norm_quad(x):
    """Normalise to quadrant -pi -> pi."""
    n = deepcopy(x)
    while np.abs(n) > np.pi:
        if n > np.pi: n -= 2*np.pi
        else: n += 2*np.pi
    return n

def Bn(B, omega, chi):
    """Normalised B."""
    nchi = norm_quad(chi)
    nomega = norm_quad(omega)
    if np.abs(nchi - nomega) >= (np.pi / 2.):
        return np.inf
    return B / np.cos(nchi - omega)

def calc_unitarity_bounds(f_s, n_samples):
    """Calculate unitarity boundary for a given source flavour ratio.

    Parameters
    ----------
    f_s : list, length = 3
        f_e, f_mu and f_tau

    n_samples : int
        number of points to sample

    Returns
    ----------
    numpy ndarray measured flavour ratios of shape (n_samples, 3)

    """
    omega = np.linspace(-np.pi, np.pi, n_samples)
    chi = np.linspace(-np.pi, np.pi, n_samples)

    x = (1 - f_s[0] - 2*f_s[1]) * np.sin(omega)
    y = (1 - 2*f_s[0] - f_s[1]) * np.cos(omega)
    z = (f_s[1] - f_s[0]) * (np.cos(omega) - np.sin(omega))

    B = B_v(x, y, z)
    if np.any(~np.isfinite(B)):
        print 'B', B
        raise AssertionError('inf elements found!')

    eta = np.full_like(chi, np.inf, dtype=np.float)
    for i_chi in xrange(n_samples):
        nB = []
        for i_ome in xrange(n_samples):
            nB.append(Bn(B[i_ome], omega[i_ome], chi[i_chi]))
        eta[i_chi] = np.min(nB)
    if np.any(~np.isfinite(eta)):
        print 'eta', eta
        raise AssertionError('inf elements found!')

    df_em = eta * np.cos(chi)
    df_um = eta * np.sin(chi)

    af_m = np.dstack([
        df_em + f_s[0],
        df_um + f_s[1],
        1 - ((df_em + f_s[0]) + (df_um + f_s[1]))
    ])[0]
    return af_m
