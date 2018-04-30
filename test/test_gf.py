import GolemFitPy as gf

def steering_params():
    steering_categ = 'p2_0'
    params = gf.SteeringParams(gf.sampleTag.HESE)
    params.fastmode = True
    params.quiet = False
    params.simToLoad= steering_categ.lower()
    return params

def set_up_as(fitter, params):
    print 'Injecting the model', params
    asimov_params = gf.FitParameters(gf.sampleTag.HESE)
    for x in params.iterkeys():
        asimov_params.__setattr__(x, float(params[x]))
    fitter.SetupAsimov(asimov_params)

def setup_fitter(asimov_paramset):
    datapaths = gf.DataPaths()
    sparams = steering_params()
    npp = gf.NewPhysicsParams()
    fitter = gf.GolemFit(datapaths, sparams, npp)
    set_up_as(fitter, asimov_paramset)
    return fitter

def get_llh(fitter, params):
    fitparams = gf.FitParameters(gf.sampleTag.HESE)
    for x in params.iterkeys():
        fitparams.__setattr__(x, float(params[x]))
    llh = -fitter.EvalLLH(fitparams)
    return llh

asimov_paramset = {'astroFlavorAngle1': 4/9., 'astroFlavorAngle2': 0.}
print 'injecting', asimov_paramset
fitter = setup_fitter(asimov_paramset)

test_paramset = {'astroFlavorAngle1': 0.36, 'astroFlavorAngle2': -0.57}
print 'testing', test_paramset
print 'llh', get_llh(fitter, test_paramset)

test_paramset = {'astroFlavorAngle1': 0.385224559219, 'astroFlavorAngle2': -0.157617854374}
print 'testing', test_paramset
print 'llh', get_llh(fitter, test_paramset)

test_paramset = {'astroFlavorAngle1': 0.415578500878, 'astroFlavorAngle2': -0.0196186993217}
print 'testing', test_paramset
print 'llh', get_llh(fitter, test_paramset)

test_paramset = {'astroFlavorAngle1': 4/9., 'astroFlavorAngle2': 0}
print 'testing', test_paramset
print 'llh', get_llh(fitter, test_paramset)

test_paramset = {'astroFlavorAngle1': 4/9, 'astroFlavorAngle2': 0}
print 'testing', test_paramset
print 'llh', get_llh(fitter, test_paramset)


import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

shape = (10, 100)
angle1_binning = np.linspace(0, 1, shape[0])
angle2_binning = np.linspace(-1, 1, shape[1])

for an1 in angle1_binning:
    llhs = []
    for an2 in angle2_binning:
        test_paramset = {'astroFlavorAngle1': an1, 'astroFlavorAngle2': an2}
        llhs.append(get_llh(fitter, test_paramset))
    plt.plot(angle2_binning, llhs, label='astroFlavorAngle1 = {0}'.format(an1))
plt.xlabel('astroFlavorAngle2')
plt.ylabel('LLH')
plt.legend()
plt.savefig('llh_profile_fastmode.png'.format(an1))
