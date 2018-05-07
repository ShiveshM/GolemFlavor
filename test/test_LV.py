#!/usr/bin/env python

from __future__ import absolute_import, division

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc

import GolemFitPy as gf

rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern'], 'size':18})

dp = gf.DataPaths()
steer = gf.SteeringParams()

fig = plt.figure(figsize=[12, 8])
ax = fig.add_subplot(111)
ax.set_xscale('log')
ax.set_yscale('log')

npp = gf.NewPhysicsParams()
npp.type = gf.NewPhysicsType.None
# npp.n_lv = 1
# npp.lambda_1 = 1.e100
# npp.lambda_2 = 1.e100

golem = gf.GolemFit(dp, steer, npp)

binning = golem.GetEnergyBinsMC()
ax.set_xlim(binning[0], binning[-1])
# ax.set_ylim(binning[0], binning[-1])

fit_params = gf.FitParameters(gf.sampleTag.HESE)
golem.SetupAsimov(fit_params)

exp = np.sum(golem.GetExpectation(fit_params), axis=(0, 1, 2, 3))
ax.step(binning, np.concatenate([[exp[0]], exp]), alpha=1,
        drawstyle='steps-pre', label='NULL', linestyle='--')

print 'NULL min_llh', golem.MinLLH().likelihood
print 'NULL expectation', exp
print

npp = gf.NewPhysicsParams()
npp.type = gf.NewPhysicsType.LorentzViolation
npp.n_lv = 1
npp.lambda_1 = 1.e20
npp.lambda_2 = 1.e20

golem.SetNewPhysicsParams(npp)

exp = np.sum(golem.GetExpectation(fit_params), axis=(0, 1, 2, 3))
ax.step(binning, np.concatenate([[exp[0]], exp]), alpha=1,
        drawstyle='steps-pre', label='1e-20 LV', linestyle='--')

print '1e20 LV min_llh', golem.MinLLH().likelihood
print '1e20 LV expectation', exp

npp = gf.NewPhysicsParams()
npp.type = gf.NewPhysicsType.LorentzViolation
npp.n_lv = 1
npp.lambda_1 = 1.e10
npp.lambda_2 = 1.e10

golem.SetNewPhysicsParams(npp)

exp = np.sum(golem.GetExpectation(fit_params), axis=(0, 1, 2, 3))
ax.step(binning, np.concatenate([[exp[0]], exp]), alpha=1,
        drawstyle='steps-pre', label='1e-10 LV', linestyle='--')

print '1e10 LV min_llh', golem.MinLLH().likelihood
print '1e10 LV expectation', exp

npp = gf.NewPhysicsParams()
npp.type = gf.NewPhysicsType.LorentzViolation
npp.n_lv = 1
npp.lambda_1 = 1.e-20
npp.lambda_2 = 1.e-20

golem.SetNewPhysicsParams(npp)

ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
ax.set_xlabel(r'Deposited energy / GeV')
ax.set_ylabel(r'Events')
for xmaj in ax.xaxis.get_majorticklocs():
    ax.axvline(x=xmaj, ls=':', color='gray', alpha=0.7, linewidth=1)
for ymaj in ax.yaxis.get_majorticklocs():
    ax.axhline(y=ymaj, ls=':', color='gray', alpha=0.7, linewidth=1)

legend = ax.legend(prop=dict(size=12))
fig.savefig('test.png', bbox_inches='tight', dpi=250)

