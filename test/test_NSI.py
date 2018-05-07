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
npp = gf.NewPhysicsParams()

steer.quiet = False
steer.fastmode = True

golem = gf.GolemFit(dp, steer, npp)

fit_params = gf.FitParameters(gf.sampleTag.HESE)
golem.SetupAsimov(fit_params)

fig = plt.figure(figsize=[6, 5])
ax = fig.add_subplot(111)
ax.set_xscale('log')
ax.set_yscale('log')

binning = golem.GetEnergyBinsMC()
ax.set_xlim(binning[0], binning[-1])
# ax.set_ylim(binning[0], binning[-1])

print 'NULL min_llh', golem.MinLLH().likelihood

# exp = np.sum(golem.GetExpectation(fit_params), axis=(0, 1, 2, 3))
# ax.step(binning, np.concatenate([[exp[0]], exp]), alpha=1,
#         drawstyle='steps-pre', label='NULL', linestyle='--')

# print 'NULL expectation', exp
print

npp.type = gf.NewPhysicsType.NonStandardInteraction
npp.epsilon_mutau = 0.1
golem.SetNewPhysicsParams(npp)

print '0.1 mutau min_llh', golem.MinLLH().likelihood

# exp = np.sum(golem.GetExpectation(fit_params), axis=(0, 1, 2, 3))
# ax.step(binning, np.concatenate([[exp[0]], exp]), alpha=1,
#         drawstyle='steps-pre', label='0.1 mutau', linestyle='--')

# print '0.1 mutau expectation', exp
print

np.epsilon_mutau = 0.2
golem.SetNewPhysicsParams(npp)

print '0.2 mutau min_llh', golem.MinLLH().likelihood

# exp = np.sum(golem.GetExpectation(fit_params), axis=(0, 1, 2, 3))
# ax.step(binning, np.concatenate([[exp[0]], exp]), alpha=1,
#         drawstyle='steps-pre', label='0.2 mutau', linestyle='--')

# print '0.2 mutau expectation', exp
print

np.epsilon_mutau = 0.3
golem.SetNewPhysicsParams(npp)

print '0.3 mutau min_llh', golem.MinLLH().likelihood

# exp = np.sum(golem.GetExpectation(fit_params), axis=(0, 1, 2, 3))
# ax.step(binning, np.concatenate([[exp[0]], exp]), alpha=1,
#         drawstyle='steps-pre', label='0.3 mutau', linestyle='--')

# print '0.3 mutau expectation', exp
print

np.epsilon_mutau = 0.4
golem.SetNewPhysicsParams(npp)

print '0.4 mutau min_llh', golem.MinLLH().likelihood

# exp = np.sum(golem.GetExpectation(fit_params), axis=(0, 1, 2, 3))
# ax.step(binning, np.concatenate([[exp[0]], exp]), alpha=1,
#         drawstyle='steps-pre', label='0.4 mutau', linestyle='--')

# print '0.4 mutau expectation', exp
print

ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
ax.set_xlabel(r'Deposited energy / GeV')
ax.set_ylabel(r'Events')
for xmaj in ax.xaxis.get_majorticklocs():
    ax.axvline(x=xmaj, ls=':', color='gray', alpha=0.7, linewidth=1)
for ymaj in ax.yaxis.get_majorticklocs():
    ax.axhline(y=ymaj, ls=':', color='gray', alpha=0.7, linewidth=1)

legend = ax.legend(prop=dict(size=12))
fig.savefig('test_NSI.png', bbox_inches='tight', dpi=250)

