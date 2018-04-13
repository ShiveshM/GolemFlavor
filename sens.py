#! /usr/bin/env python
# author : S. Mandalia
#          s.p.mandalia@qmul.ac.uk
#
# date   : April 10, 2018

"""
HESE BSM flavour ratio analysis script
"""

from __future__ import absolute_import, division

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import rc

import GolemFitPy as gf

import fr
from utils import gf as gf_utils
from utils import likelihood as llh_utils
from utils import misc as misc_utils
from utils.enums import Likelihood, ParamTag
from utils.plot import plot_BSM_angles_limit

rc('text', usetex=False)
rc('font', **{'family':'serif', 'serif':['Computer Modern'], 'size':18})


RUN = False


z = 0+1E-6
scenarios = [
    [np.sin(np.pi/2.)**2, z, z, z],
    [z, np.cos(np.pi/2.)**4, z, z],
    [z, z, np.sin(np.pi/2.)**2, z],
]
xticks = [r'$\mathcal{O}_{12}$', r'$\mathcal{O}_{13}$', r'$\mathcal{O}_{23}$']

def fit_flags(flag_dict):
    flags = gf.FitParametersFlag()
    for key in flag_dict.iterkeys():
        flags.__setattr__(key, flag_dict[key])
    return flags

default_flags = {
    # False means it's not fixed in minimization
    'astroFlavorAngle1'         : True,
    'astroFlavorAngle2'         : True,
    # 'astroENorm'                : True,
    # 'astroMuNorm'               : True,
    # 'astroTauNorm'              : True,
    'convNorm'                  : False,
    'promptNorm'                : False,
    'muonNorm'                  : False,
    'astroNorm'                 : False,
    'astroParticleBalance'      : True,
    'astroDeltaGamma'           : False,
    'cutoffEnergy'              : True,
    'CRDeltaGamma'              : False,
    'piKRatio'                  : False,
    'NeutrinoAntineutrinoRatio' : True,
    'darkNorm'                  : True,
    'domEfficiency'             : True,
    'holeiceForward'            : True,
    'anisotropyScale'           : True,
    'astroNormSec'              : True,
    'astroDeltaGammaSec'        : True
}


def main():
    args = fr.parse_args()
    args.likelihood = Likelihood.GF_FREQ
    fr.process_args(args)
    misc_utils.print_args(args)

    bins = 10

    asimov_paramset, mcmc_paramset = fr.get_paramsets(args)
    outfile = misc_utils.gen_outfile_name(args)
    print '== {0:<25} = {1}'.format('outfile', outfile)

    asimov_paramset = asimov_paramset.from_tag(ParamTag.BESTFIT)
    mcmc_paramset = mcmc_paramset.from_tag(ParamTag.NUISANCE, invert=True)

    sc_range = mcmc_paramset.from_tag(ParamTag.SCALE)[0].ranges
    scan_scales = np.linspace(sc_range[0], sc_range[1], bins+1)
    print 'scan_scales', scan_scales

    outfile = './sens'
    if RUN:
        datapaths = gf.DataPaths()
        sparams = gf_utils.steering_params(args)
        npp = gf.NewPhysicsParams()
        fitter = gf.GolemFit(datapaths, sparams, npp)
        fitter.SetFitParametersFlag(fit_flags(default_flags))
        gf_utils.set_up_as(fitter, asimov_paramset)

        data = np.zeros((len(scenarios), bins+1, 2))
        mm_angles = mcmc_paramset.from_tag(ParamTag.MMANGLES)
        for idx, scen in enumerate(scenarios):
            arr = []
            scales = []
            llhs = []
            for yidx, an in enumerate(mm_angles):
                an.value = scen[yidx]
            for sc in scan_scales:
                theta = scen + [sc]
                print 'theta', theta
                llh = llh_utils.triangle_llh(
                    theta=theta, args=args, asimov_paramset=asimov_paramset,
                    mcmc_paramset=mcmc_paramset, fitter=fitter
                )
                print 'llh', llh
                scales.append(sc)
                llhs.append(llh)

            for i, d in enumerate(llhs):
                data[idx][i][0] = scales[i]
                data[idx][i][1] = d

        np.save(outfile + '.npy', data)

    plot_BSM_angles_limit(
        infile=outfile+'.npy',
        outfile=outfile,
        xticks=xticks,
        outformat=['png'],
        args=args
    )


main.__doc__ = __doc__


if __name__ == '__main__':
    main()

