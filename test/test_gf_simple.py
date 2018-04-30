import GolemFitPy as gf

dp = gf.DataPaths()
npp = gf.NewPhysicsParams()
sp = gf.SteeringParams(gf.sampleTag.HESE)

sp.quiet = False
# sp.fastmode = True

golem = gf.GolemFit(dp, sp, npp)

fp = gf.FitParameters(gf.sampleTag.HESE)
fp.astroFlavorAngle1 = 4./9.
fp.astroFlavorAngle2 = 0

golem.SetupAsimov(fp)

fp_sh = gf.FitParameters(gf.sampleTag.HESE)
fp_sh.astroFlavorAngle1 = 0.36
fp_sh.astroFlavorAngle2 = -0.57

print 'Eval fp = {0}'.format(golem.EvalLLH(fp))
print 'Eval fp_sh = {0}'.format(golem.EvalLLH(fp_sh))
