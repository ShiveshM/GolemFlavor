import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import GolemFitPy as gf

FASTMODE = True
PARAMETERS = [
    # 'astroFlavorAngle1', 'astroFlavorAngle2',
    'convNorm',
    # 'promptNorm', 'muonNorm', 'astroNorm'
]
DEFAULTS = [
    # 4/9., 0., 1., 0., 1., 6.9
    1.
]
RANGES = [
    # (0, 1), (-1, 1), (0.01, 10), (0., 30), (0.01, 10), (0.01, 30)
    (0.01, 10)
]
BINS = 50

def steering_params():
    steering_categ = 'p2_0'
    params = gf.SteeringParams(gf.sampleTag.HESE)
    if FASTMODE:
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
    priors = gf.Priors()
    priors.convNormWidth = 9e9
    fitter.SetFitPriors(priors)

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

for ip, param in enumerate(PARAMETERS):
    asimov_paramset = {param: DEFAULTS[ip]}
    print 'injecting', asimov_paramset
    fitter = setup_fitter(asimov_paramset)
    binning = np.linspace(RANGES[ip][0], RANGES[ip][1], BINS)
    llhs = []
    for b in binning:
        test_paramset = {param: b}
        print 'testing', test_paramset
        llh = get_llh(fitter, test_paramset)
        print 'llh', llh
        llhs.append(llh)
    plt.plot(binning, llhs)
    plt.axvline(x=DEFAULTS[ip])
    plt.xlabel(param)
    plt.ylabel('LLH')
    outfile = 'llh_profile_noprior_'
    if FASTMODE:
        plt.savefig(outfile + 'fastmode_{0}.png'.format(param))
    else:
        plt.savefig(outfile + '{0}.png'.format(param))
    plt.clf()
