#! /usr/bin/env python
"""
From a gaussian likelihood run an MCMC scan to find the posteriors
"""

from __future__ import absolute_import, division

import os, sys
import errno

import argparse
import multiprocessing
from copy import deepcopy

import numpy as np
from scipy import linalg

import emcee
import tqdm

import GolemFitPy as gf
import chainer_plot


RUN_SCAN = True

GOLEMFIT = None
FIX_MIXING = False
FIX_SFR = True
SOURCE_FR = [1, 2, 0]
FIX_SCALE = True
NO_BSM = False

DIMENSION = 3
ENERGY = 1000000 # GeV
MEASURED_FR = [1, 1, 1]
SIGMA = 0.001
SCALE_RATIO = 100.
MASS_EIGENVALUES = [7.40E-23, 2.515E-21]
# MASS_EIGENVALUES = [7.40E-100, 2.515E-100]
FLAT = False

CHANGE_MASS = False
if CHANGE_MASS:
    SCALE = 1E-27
else:
    SCALE = np.power(10, np.round(np.log10(MASS_EIGENVALUES[1]/ENERGY)) - np.log10(ENERGY**(DIMENSION-3)))
SCALE2_BOUNDS = (SCALE*1E-10, SCALE*1E10)


def set_up_gf():
    steering_params = gf.SteeringParams()
    datapaths = gf.DataPaths()
    npp = gf.NewPhysicsParams()

    steering_params.logEbinWidth=(7.0-4.78)/20.
    steering_params.cosThbinWidth=0.2
    steering_params.correct_prompt_knee=True
    steering_params.fastmode=True
    # TODO(shivesh): figure out right number for missid
    steering_params.track_to_shower_missid = 0.3
    golemfitsourcepath = os.environ['GOLEMSOURCEPATH'] + '/GolemFit'
    datapaths.compact_file_path = golemfitsourcepath + '/compact_data/'  # only used for FastMode; not currently working
    datapaths.data_path = golemfitsourcepath + '/data/'
    datapaths.mc_path = golemfitsourcepath + '/monte_carlo/'

    golemfit = gf.GolemFit(datapaths, steering_params, npp)
    return golemfit


def set_up_as(golemfit, parms):
    asimov_parameters = gf.FitParameters()
    asimov_parameters.convNorm = parms['convNorm']
    asimov_parameters.promptNorm = parms['promptNorm']
    asimov_parameters.muonNorm = parms['muonNorm']
    asimov_parameters.astroNorm = parms['astroNorm']
    asimov_parameters.astroDeltaGamma = parms['astroDeltaGamma']
    asimov_parameters.astroENorm = parms['astroENorm']
    asimov_parameters.astroMuNorm = parms['astroMuNorm']
    asimov_parameters.astroTauNorm = parms['astroTauNorm']
    golemfit.SetupAsimov(asimov_parameters)


def get_llh(golemfit, parms):
    fitparams = gf.FitParameters()
    fitparams.convNorm = parms['convNorm']
    fitparams.promptNorm = parms['promptNorm']
    fitparams.muonNorm = parms['muonNorm']
    fitparams.astroNorm = parms['astroNorm']
    fitparams.astroDeltaGamma = parms['astroDeltaGamma']
    fitparams.astroENorm = parms['astroENorm']
    fitparams.astroMuNorm = parms['astroMuNorm']
    fitparams.astroTauNorm = parms['astroTauNorm']
    return golemfit.EvalLLH(fitparams)


def test_uni(x):
    """Test the unitarity of a matrix."""
    # print 'Unitarity test:\n{0}'.format(abs(np.dot(x, x.conj().T)))
    return abs(np.dot(x, x.conj().T))


def angles_to_fr(angles):
    sphi4, c2psi = angles

    psi = (0.5)*np.arccos(c2psi)

    sphi2 = np.sqrt(sphi4)
    cphi2 = 1. - sphi2
    spsi2 = np.sin(psi)**2
    cspi2 = 1. - spsi2

    x = sphi2*cspi2
    y = sphi2*spsi2
    z = cphi2
    return x, y, z


def angles_to_u(angles):
    s12_2, c13_4, s23_2, dcp = angles
    dcp = np.complex128(dcp)

    c12_2 = 1. - s12_2
    c13_2 = np.sqrt(c13_4)
    s13_2 = 1. - c13_2
    c23_2 = 1. - s23_2

    t12 = np.arcsin(np.sqrt(s12_2))
    t13 = np.arccos(np.sqrt(c13_2))
    t23 = np.arcsin(np.sqrt(s23_2))

    c12 = np.cos(t12)
    s12 = np.sin(t12)
    c13 = np.cos(t13)
    s13 = np.sin(t13)
    c23 = np.cos(t23)
    s23 = np.sin(t23)

    p1 = np.array([[1   , 0   , 0]                   , [0    , c23 , s23] , [0                   , -s23 , c23]] , dtype=np.complex128)
    p2 = np.array([[c13 , 0   , s13*np.exp(-1j*dcp)] , [0    , 1   , 0]   , [-s13*np.exp(1j*dcp) , 0    , c13]] , dtype=np.complex128)
    p3 = np.array([[c12 , s12 , 0]                   , [-s12 , c12 , 0]   , [0                   , 0    , 1]]   , dtype=np.complex128)

    u = np.dot(np.dot(p1, p2), p3)
    return u


def cardano_eqn(ham):
    """Diagonalise the effective Hamiltonian 3x3 matrix into the form
    h_{eff} = UE_{eff}U^{dagger} using the procedure in PRD91, 052003 (2015)
    """
    a = -np.trace(ham)
    b = (0.5) * ((np.trace(ham))**2 - np.trace(np.dot(ham, ham)))
    c = -linalg.det(ham)

    Q = (1/9.) * (a**2 - 3*b)
    R = (1/54.) * (2*a**3 - 9*a*b + 27*c)
    theta = np.arccos(R / np.sqrt(Q**3))

    E1 = -2 * np.sqrt(Q) * np.cos(theta/3.) - (1/3.)*a
    E2 = -2 * np.sqrt(Q) * np.cos((theta - 2.*np.pi)/3.) - (1/3.)*a
    E3 = -2 * np.sqrt(Q) * np.cos((theta + 2.*np.pi)/3.) - (1/3.)*a

    A1 = ham[1][2] * (ham[0][0] - E1) - ham[1][0]*ham[0][2]
    A2 = ham[1][2] * (ham[0][0] - E2) - ham[1][0]*ham[0][2]
    A3 = ham[1][2] * (ham[0][0] - E3) - ham[1][0]*ham[0][2]

    B1 = ham[2][0] * (ham[1][1] - E1) - ham[2][1]*ham[1][0]
    B2 = ham[2][0] * (ham[1][1] - E2) - ham[2][1]*ham[1][0]
    B3 = ham[2][0] * (ham[1][1] - E3) - ham[2][1]*ham[1][0]

    C1 = ham[1][0] * (ham[2][2] - E1) - ham[1][2]*ham[2][0]
    C2 = ham[1][0] * (ham[2][2] - E2) - ham[1][2]*ham[2][0]
    C3 = ham[1][0] * (ham[2][2] - E3) - ham[1][2]*ham[2][0]

    N1 = np.sqrt(abs(A1*B1)**2 + abs(A1*C1)**2 + abs(B1*C1)**2)
    N2 = np.sqrt(abs(A2*B2)**2 + abs(A2*C2)**2 + abs(B2*C2)**2)
    N3 = np.sqrt(abs(A3*B3)**2 + abs(A3*C3)**2 + abs(B3*C3)**2)

    mm = np.array([
        [np.conjugate(B1)*C1 / N1, np.conjugate(B2)*C2 / N2, np.conjugate(B3)*C3 / N3],
        [A1*C1 / N1, A2*C2 / N2, A3*C3 / N3],
        [A1*B1 / N1, A2*B2 / N2, A3*B3 / N3]
    ])
    return mm


# s_12^2, c_13^4, s_23^2, dcp
NUFIT_U = angles_to_u((0.307, (1-0.02195)**2, 0.565, 3.97935))


def params_to_BSMu(theta):
    if FIX_MIXING:
        s12_2, c13_4, s23_2, dcp, sc2 = 0.5, 1.0-1E-6, 0.5, 0., theta
    elif FIX_SCALE:
        s12_2, c13_4, s23_2, dcp = theta
        sc2 = np.log10(SCALE)
    else:
        s12_2, c13_4, s23_2, dcp, sc2 = theta
    sc2 = np.power(10., sc2)
    sc1 = sc2 / SCALE_RATIO

    mass_matrix = np.array(
        [[0, 0, 0], [0, MASS_EIGENVALUES[0], 0], [0, 0, MASS_EIGENVALUES[1]]]
    )
    sm_ham = (1./(2*ENERGY))*np.dot(NUFIT_U, np.dot(mass_matrix, NUFIT_U.conj().T))
    if NO_BSM:
        eg_vector = cardano_eqn(sm_ham)
        tu = test_uni(eg_vector)
        if not abs(np.trace(tu) - 3.) < 1e-5 or not abs(np.sum(tu) - 3.) < 1e-5:
            raise AssertionError('Matrix is not unitary!\neg_vector\n{0}\ntest u\n{1}'.format(eg_vector, tu))
        return eg_vector

    new_physics_u = angles_to_u((s12_2, c13_4, s23_2, dcp))
    scale_matrix = np.array(
        [[0, 0, 0], [0, sc1, 0], [0, 0, sc2]]
    )
    bsm_term = (ENERGY**(DIMENSION-3)) * np.dot(new_physics_u, np.dot(scale_matrix, new_physics_u.conj().T))

    bsm_ham = sm_ham + bsm_term

    eg_vector = cardano_eqn(bsm_ham)
    tu = test_uni(eg_vector)
    if not abs(np.trace(tu) - 3.) < 1e-5 or not abs(np.sum(tu) - 3.) < 1e-5:
	raise AssertionError('Matrix is not unitary!\neg_vector\n{0}\ntest u\n{1}'.format(eg_vector, tu))
    return eg_vector


def u_to_fr(initial_fr, matrix):
    """Compute the observed flavour ratio assuming decoherence."""
    # TODO(shivesh): energy dependence
    composition = np.einsum(
        'ai, bi, a -> b', abs(matrix)**2, abs(matrix)**2, initial_fr
    )
    ratio = composition / np.sum(initial_fr)
    return ratio


def triangle_llh(theta):
    """-Log likelihood function for a given theta."""
    convNorm, promptNorm, muonNorm, astroNorm, astroDeltaGamma = theta[:5]
    theta = deepcopy(theta[5:])

    if FIX_SFR:
        fr1, fr2, fr3 = SOURCE_FR
        u = params_to_BSMu(theta)
    else:
        fr1, fr2, fr3 = angles_to_fr(theta[-2:])
        u = params_to_BSMu(theta[:-2])

    fr = u_to_fr((fr1, fr2, fr3), u)

    fit_values = {
        'convNorm': convNorm,
        'promptNorm': promptNorm,
        'muonNorm': muonNorm,
        'astroNorm': astroNorm,
        'astroDeltaGamma': astroDeltaGamma,
        'astroENorm': fr[0],
        'astroMuNorm': fr[1],
        'astroTauNorm': fr[2]
    }
    if FLAT:
        return 10.
    else:
        return get_llh(GOLEMFIT, fit_values)


def lnprior(theta):
    """Priors on theta."""
    convNorm, promptNorm, muonNorm, astroNorm, astroDeltaGamma = theta[:5]
    theta = deepcopy(theta[5:])

    # Nuisance parameter bounds
    if 0. < convNorm < 5. and 0. < promptNorm < 5. and 0. < muonNorm < 5. and \
       0. < astroNorm < 5. and 0. < astroDeltaGamma < 5.:
        pass
    else: return -np.inf

    if FIX_SFR:
        if FIX_MIXING:
            s12_2, sc2 = theta
        elif FIX_SCALE:
            s12_2, c13_4, s23_2, dcp = theta
        else:
            s12_2, c13_4, s23_2, dcp, sc2 = theta
    else:
        if FIX_MIXING:
            sc2, sphi4, c2psi = theta
        elif FIX_SCALE:
            s12_2, c13_4, s23_2, dcp, sphi4, c2psi = theta
        else:
            s12_2, c13_4, s23_2, dcp, sc2, sphi4, c2psi = theta
    if not FIX_SCALE:
        sc2 = np.power(10., sc2)

    if not FIX_SFR:
        # Flavour ratio bounds
        if 0. <= sphi4 <= 1.0 and -1.0 <= c2psi <= 1.0:
            pass
        else: return -np.inf

    # Mixing angle bounds
    if not FIX_MIXING:
        if 0. <= s12_2 <= 1. and 0. <= c13_4 <= 1. and 0. <= s23_2 <= 1. \
           and 0. <= dcp <= 2*np.pi:
            pass
        else: return -np.inf

    # Scale bounds
    if not FIX_SCALE:
        if SCALE2_BOUNDS[0] <= sc2 <= SCALE2_BOUNDS[1]:
            pass
        else: return -np.inf

    return 0.


def lnprob(theta):
    """Prob function for mcmc."""
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + triangle_llh(theta)


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--measured-ratio', type=int, nargs=3, default=[1, 1, 1],
        help='Set the central value for the measured flavour ratio at IceCube'
    )
    parser.add_argument(
        '--sigma-ratio', type=float, default=0.01,
        help='Set the 1 sigma for the measured flavour ratio'
    )
    parser.add_argument(
        '--fix-source-ratio', type=str, default='False',
        help='Fix the source flavour ratio'
    )
    parser.add_argument(
        '--source-ratio', type=int, nargs=3, default=[2, 1, 0],
        help='Set the source flavour ratio for the case when you want to fix it'
    )
    parser.add_argument(
        '--no-bsm', type=str, default='False',
        help='Turn off BSM terms'
    )
    parser.add_argument(
        '--fix-mixing', type=str, default='False',
        help='Fix all mixing parameters except one'
    )
    parser.add_argument(
        '--fix-scale', type=str, default='False',
        help='Fix the new physics scale'
    )
    parser.add_argument(
        '--scale', type=float, default=1e-30,
        help='Set the new physics scale'
    )
    parser.add_argument(
        '--dimension', type=int, default=3,
        help='Set the new physics dimension to consider'
    )
    parser.add_argument(
        '--energy', type=float, default=1000,
        help='Set the energy scale'
    )
    parser.add_argument(
        '--flat-llh', type=str, default='False',
        help='Run with a flat likelihood'
    )
    parser.add_argument(
        '--burnin', type=int, default=100,
        help='Amount to burnin'
    )
    parser.add_argument(
        '--nwalkers', type=int, default=100,
        help='Number of walkers'
    )
    parser.add_argument(
        '--nsteps', type=int, default=2000,
        help='Number of steps to run'
    )
    parser.add_argument(
        '--seed', type=int, default=99,
        help='Set the random seed value'
    )
    parser.add_argument(
        '--outfile', type=str, default='./untitled',
        help='Path to output chains'
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    np.random.seed(args.seed)

    global DIMENSION
    global ENERGY
    global MEASURED_FR
    global SIGMA
    global FLAT
    global FIX_SFR
    global SOURCE_FR
    global FIX_MIXING
    global FIX_SCALE
    global SCALE
    global SCALE2_BOUNDS
    global NO_BSM

    DIMENSION = args.dimension
    ENERGY = args.energy

    MEASURED_FR = np.array(args.measured_ratio) / float(np.sum(args.measured_ratio))
    SIGMA = args.sigma_ratio

    if args.flat_llh.lower() == 'true':
        FLAT = True
    elif args.flat_llh.lower() == 'false':
        FLAT = False
    else:
        raise ValueError

    if args.fix_source_ratio.lower() == 'true':
        FIX_SFR = True
    elif args.fix_source_ratio.lower() == 'false':
        FIX_SFR = False
    else:
        raise ValueError

    if args.fix_mixing.lower() == 'true':
        FIX_MIXING = True
    elif args.fix_mixing.lower() == 'false':
        FIX_MIXING = False
    else:
        raise ValueError

    if args.fix_scale.lower() == 'true':
        FIX_SCALE = True
    elif args.fix_scale.lower() == 'false':
        FIX_SCALE = False
    else:
        raise ValueError

    if args.no_bsm.lower() == 'true':
        NO_BSM = True
    elif args.no_bsm.lower() == 'false':
        NO_BSM = False
    else:
        raise ValueError

    if FIX_SFR:
        SOURCE_FR = np.array(args.source_ratio) / float(np.sum(args.source_ratio))

    if FIX_SCALE:
	SCALE = args.scale
    else:
        if CHANGE_MASS:
            SCALE = 1E-27
        else:
            SCALE = np.power(10, np.round(np.log10(MASS_EIGENVALUES[1]/ENERGY)) - np.log10(ENERGY**(DIMENSION-3)))

    if FIX_MIXING and FIX_SFR:
        raise NotImplementedError('Fixed mixing and sfr not implemented')
    if FIX_MIXING and FIX_SCALE:
        raise NotImplementedError('Fixed mixing and scale not implemented')

    SCALE2_BOUNDS = (SCALE*1E-10, SCALE*1E10)

    print 'RUN_SCAN = {0}'.format(RUN_SCAN)
    print 'MEASURED_FR = {0}'.format(MEASURED_FR)
    print 'SIGMA = {0}'.format(SIGMA)
    print 'FLAT = {0}'.format(FLAT)
    print 'NO_BSM = {0}'.format(NO_BSM)
    print 'ENERGY = {0}'.format(ENERGY)
    print 'DIMENSION = {0}'.format(DIMENSION)
    print 'SCALE2_BOUNDS = {0}'.format(SCALE2_BOUNDS)
    print 'FIX_SFR = {0}'.format(FIX_SFR)
    if FIX_SFR:
        print 'SOURCE_FR = {0}'.format(SOURCE_FR)
    print 'FIX_MIXING = {0}'.format(FIX_MIXING)
    print 'FIX_SCALE = {0}'.format(FIX_SCALE)
    if FIX_SCALE:
        print 'SCALE = {0}'.format(SCALE)

    global GOLEMFIT
    if RUN_SCAN:
        GOLEMFIT = set_up_gf()

        # TODO(shivesh): make all these into args
        inject = {
            'convNorm': 1.,
            'promptNorm': 0.,
            'muonNorm': 1.,
            'astroNorm': 1.,
            'astroDeltaGamma': 2.,
            'astroENorm': MEASURED_FR[0],
            'astroMuNorm': MEASURED_FR[1],
            'astroTauNorm': MEASURED_FR[2],
        }

        print "Injecting this model: ", inject
        set_up_as(GOLEMFIT, inject)

    if FIX_SFR:
        if FIX_MIXING:
            ndim = 7
        elif FIX_SCALE:
            ndim = 9
        else:
            ndim = 10
    else:
        if FIX_MIXING:
            ndim = 8
        elif FIX_SCALE:
            ndim = 11
        else:
            ndim = 12
    nwalkers = args.nwalkers
    ntemps = 1
    burnin = args.burnin
    betas = np.array([1e0, 1e-1, 1e-2, 1e-3, 1e-4])
    s2 = np.average(np.log10(SCALE2_BOUNDS))
    if FIX_SFR:
        if FIX_MIXING:
            p0_base = [1., 0., 1., 1., 2., 0.5, s2]
            p0_std = [0.3, 0.05, 0.1, 0.1, 0.1, 0.2, 3]
        elif FIX_SCALE:
            p0_base = [1., 0., 1., 1., 2., 0.5, 0.5, 0.5, np.pi]
            p0_std = [0.3, 0.05, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2]
        else:
            p0_base = [1., 0., 1., 1., 2., 0.5, 0.5, 0.5, np.pi, s2]
            p0_std = [0.3, 0.05, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 3]
    else:
        if FIX_MIXING:
            p0_base = [1., 0., 1., 1., 2., s2, 0.5, 0.0]
            p0_std = [0.3, 0.05, 0.1, 0.1, 0.1, 3, 0.2, 0.2]
        elif FIX_SCALE:
            p0_base = [1., 0., 1., 1., 2., 0.5, 0.5, 0.5, np.pi, 0.5, 0.0]
            p0_std = [0.3, 0.05, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
        else:
            p0_base = [1., 0., 1., 1., 2., 0.5, 0.5, 0.5, np.pi, s2, 0.5, 0.0]
            p0_std = [0.3, 0.05, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 3, 0.2, 0.2]

    print 'p0_base', p0_base
    print 'p0_std', p0_std
    p0 = np.random.normal(p0_base, p0_std, size=[ntemps, nwalkers, ndim])
    print map(lnprior, p0[0])

    if RUN_SCAN:
        threads = multiprocessing.cpu_count()
        # threads = 1
        sampler = emcee.PTSampler(
            ntemps, nwalkers, ndim, triangle_llh, lnprior, threads=threads
        )

    if RUN_SCAN:
        print "Running burn-in"
        for result in tqdm.tqdm(sampler.sample(p0, iterations=burnin), total=burnin):
            pos, prob, state = result
        sampler.reset()
        print "Finished burn-in"

        nsteps = args.nsteps

        print "Running"
        for _ in tqdm.tqdm(sampler.sample(pos, iterations=nsteps), total=nsteps):
            pass
        print "Finished"

    if FIX_SFR:
        if FIX_MIXING:
            outfile = args.outfile+'_{0:03d}_{1:03d}_{2:03d}_{3:04d}_sfr_{4:03d}_{5:03d}_{6:03d}_DIM{7}_fix_mixing'.format(
                int(MEASURED_FR[0]*100), int(MEASURED_FR[1]*100), int(MEASURED_FR[2]*100), int(SIGMA*1000),
                int(SOURCE_FR[0]*100), int(SOURCE_FR[1]*100), int(SOURCE_FR[2]*100), DIMENSION
            )
        elif FIX_SCALE:
            outfile = args.outfile+'_{0:03d}_{1:03d}_{2:03d}_{3:04d}_sfr_{4:03d}_{5:03d}_{6:03d}_DIM{7}_fixed_scale_{8}'.format(
                int(MEASURED_FR[0]*100), int(MEASURED_FR[1]*100), int(MEASURED_FR[2]*100), int(SIGMA*1000),
                int(SOURCE_FR[0]*100), int(SOURCE_FR[1]*100), int(SOURCE_FR[2]*100), DIMENSION, SCALE
            )
        else:
            outfile = args.outfile+'_{0:03d}_{1:03d}_{2:03d}_{3:04d}_sfr_{4:03d}_{5:03d}_{6:03d}_DIM{7}_single_scale'.format(
                int(MEASURED_FR[0]*100), int(MEASURED_FR[1]*100), int(MEASURED_FR[2]*100), int(SIGMA*1000),
                int(SOURCE_FR[0]*100), int(SOURCE_FR[1]*100), int(SOURCE_FR[2]*100), DIMENSION
            )
    else:
        if FIX_MIXING:
            outfile = args.outfile+'_{0:03d}_{1:03d}_{2:03d}_{3:04d}_DIM{4}_fix_mixing'.format(
                int(MEASURED_FR[0]*100), int(MEASURED_FR[1]*100), int(MEASURED_FR[2]*100),
                int(SIGMA*1000), DIMENSION
            )
        elif FIX_SCALE:
            outfile = args.outfile+'_{0:03d}_{1:03d}_{2:03d}_{3:04d}_DIM{4}_fixed_scale_{5}'.format(
                int(MEASURED_FR[0]*100), int(MEASURED_FR[1]*100), int(MEASURED_FR[2]*100),
                int(SIGMA*1000), DIMENSION, SCALE
            )
        else:
            outfile = args.outfile+'_{0:03d}_{1:03d}_{2:03d}_{3:04d}_DIM{4}'.format(
                int(MEASURED_FR[0]*100), int(MEASURED_FR[1]*100), int(MEASURED_FR[2]*100),
                int(SIGMA*1000), DIMENSION
            )

    if RUN_SCAN:
        samples = sampler.chain[0, :, :, :].reshape((-1, ndim))
        print 'acceptance fraction', sampler.acceptance_fraction
        print 'sum of acceptance fraction', np.sum(sampler.acceptance_fraction)
        print 'np.unique(samples[:,0]).shape', np.unique(samples[:,0]).shape

        try:
            print 'autocorrelation', sampler.acor
        except:
            print 'WARNING : NEED TO RUN MORE SAMPLES FOR FILE ' + outfile

    if FLAT:
	outfile += '_flat'

    if RUN_SCAN:
        try:
            os.makedirs(outfile[:-len(os.path.basename(outfile))])
        except OSError as exc:  # Python >2.5
            if exc.errno == errno.EEXIST and os.path.isdir(outfile[:-len(os.path.basename(outfile))]):
                pass
            else:
                raise
        np.save(outfile+'.npy', samples)

    print "Making triangle plots"
    chainer_plot.plot(
        infile=outfile+'.npy',
        angles=True,
        outfile=outfile[:5]+outfile[5:].replace('data', 'plots')+'_angles.pdf',
        measured_ratio=MEASURED_FR,
        sigma_ratio=SIGMA,
        fix_sfr=FIX_SFR,
        fix_mixing=FIX_MIXING,
        fix_scale=FIX_SCALE,
        source_ratio=SOURCE_FR,
        scale=SCALE,
	dimension=DIMENSION,
	energy=ENERGY,
	scale_bounds=SCALE2_BOUNDS,
    )
    # chainer_plot.plot(
    #     infile=outfile+'.npy',
    #     angles=False,
    #     outfile=outfile[:5]+outfile[5:].replace('data', 'plots')+'.pdf',
    #     measured_ratio=MEASURED_FR,
    #     sigma_ratio=SIGMA,
    #     fix_sfr=FIX_SFR,
    #     fix_mixing=FIX_MIXING,
    #     fix_scale=FIX_SCALE,
    #     source_ratio=SOURCE_FR,
    #     scale=SCALE,
    #     dimension=DIMENSION,
    #     energy=ENERGY,
    #     scale_bounds=SCALE2_BOUNDS,
    # )
    print "DONE!"


main.__doc__ = __doc__


if __name__ == '__main__':
    main()

