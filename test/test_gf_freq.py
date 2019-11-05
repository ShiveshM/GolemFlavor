import GolemFitPy as gf

FASTMODE = False

def steering_params():
    steering_categ = 'p2_0'
    params = gf.SteeringParams(gf.sampleTag.MagicTau)
    params.quiet = False
    if FASTMODE:
        params.fastmode = True
    else:
        params.fastmode = False
    params.simToLoad= steering_categ.lower()
    params.evalThreads = 4
    params.minFitEnergy = 6E4 # GeV
    params.maxFitEnergy = 1E7 # GeV
    params.load_data_from_text_file = False
    params.do_HESE_reshuffle=False
    params.use_legacy_selfveto_calculation = False
    return params

def setup_fitter():
    datapaths = gf.DataPaths()
    sparams = steering_params()
    npp = gf.NewPhysicsParams()
    fitter = gf.GolemFit(datapaths, sparams, npp)
    return fitter

def fit_flags(fitter):
    default_flags = {
        # False means it's not fixed in minimization
        'astroFlavorAngle1'         : False,
        'astroFlavorAngle2'         : False,
        'astroNorm'                 : True,
    }
    flags = gf.FitParametersFlag()
    gf_nuisance = []
    for param in default_flags.keys():
        if default_flags[param]:
            flags.__setattr__(param, True)
        else:
            print 'Setting param {0:<15} to float in the ' \
                'minimisation'.format(param)
            flags.__setattr__(param, False)
    fitter.SetFitParametersFlag(flags)

def set_up_as(fitter, params):
    print 'Injecting the model', params
    asimov_params = gf.FitParameters(gf.sampleTag.MagicTau)
    for x in params.keys():
        asimov_params.__setattr__(x, float(params[x]))
    fitter.SetFitParametersSeed(asimov_params)

def get_bf_freq(fitter):
    bf = fitter.MinLLH()
    return bf

# Setup fitter
fitter = setup_fitter()
fit_flags(fitter)

params = {'astroFlavorAngle1': 4/9., 'astroFlavorAngle2': 0.}
print
set_up_as(fitter, params)
print 'fitting...'
bf = get_bf_freq(fitter)
print 'bestfit params = astroFlavorAngle1:', bf.params.astroFlavorAngle1, \
      ', astroFlavorAngle2:', bf.params.astroFlavorAngle2
print 'bestfit llh =', -bf.likelihood
print

params = {'astroFlavorAngle1': 2/6., 'astroFlavorAngle2': 1/2.}
print
set_up_as(fitter, params)
print 'fitting...'
bf = get_bf_freq(fitter)
print 'bestfit params = astroFlavorAngle1:', bf.params.astroFlavorAngle1, \
      ', astroFlavorAngle2:', bf.params.astroFlavorAngle2
print 'bestfit llh =', -bf.likelihood
print
