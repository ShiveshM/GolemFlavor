import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import GolemFitPy as gf

FASTMODE = False

dp = gf.DataPaths()
npp = gf.NewPhysicsParams()
sp = gf.SteeringParams(gf.sampleTag.MagicTau)

sp.quiet = False
if FASTMODE:
    sp.fastmode = True
# sp.frequentist = True
sp.load_data_from_text_file = False

golem = gf.GolemFit(dp, sp, npp)

fp = gf.FitParameters(gf.sampleTag.MagicTau)
fp.astroFlavorAngle1 = 4./9.
fp.astroFlavorAngle2 = 0.

# golem.SetupAsimov(fp)
seed = 0
golem.Swallow(golem.SpitRealization(fp, seed))

fp_sh = gf.FitParameters(gf.sampleTag.MagicTau)
# fp_sh.astroFlavorAngle1 = 0.36
# fp_sh.astroFlavorAngle2 = -0.57
fp_sh.astroFlavorAngle1 = 0.
fp_sh.astroFlavorAngle2 = 1.

print 'Eval fp = {0}'.format(golem.EvalLLH(fp))

# energy_centers = golem.GetEnergyBinsMC()[:-1]+ np.diff(golem.GetEnergyBinsMC())/2.

# plt.hist(energy_centers,bins=golem.GetEnergyBinsMC(),
#          weights=np.sum(golem.GetExpectation(fp),axis=(0,1,2,3)),
#          histtype="step", lw = 2, label='injected')

# data_energy_dist = np.sum(golem.GetDataDistribution(),axis=(0,1,2,3))
# energy_centers=golem.GetEnergyBinsData()[:-1]+ np.diff(golem.GetEnergyBinsData())/2.
# plt.errorbar(energy_centers,data_energy_dist,yerr = np.sqrt(data_energy_dist),fmt='o')

print 'Eval fp_sh = {0}'.format(golem.EvalLLH(fp_sh))

# plt.hist(energy_centers,bins=golem.GetEnergyBinsMC(),
#          weights=np.sum(golem.GetExpectation(fp_sh),axis=(0,1,2,3)),
#          histtype="step", lw = 2, label='test')

# data_energy_dist = np.sum(golem.GetDataDistribution(),axis=(0,1,2,3))
# energy_centers=golem.GetEnergyBinsData()[:-1]+ np.diff(golem.GetEnergyBinsData())/2.
# plt.errorbar(energy_centers,data_energy_dist,yerr = np.sqrt(data_energy_dist),fmt='o')

# plt.loglog(nonposy="clip")
# plt.xlabel(r"Deposited energy/GeV")
# plt.ylabel(r"Events")

# outname = 'Expectation'
# if FASTMODE:
#     plt.savefig(outname + 'fastmode.png')
# else:
#     plt.savefig(outname + '.png')
