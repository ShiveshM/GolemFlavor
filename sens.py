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
import pymultinest

import fr
from utils import gf as gf_utils
from utils import likelihood as llh_utils
from utils import misc as misc_utils
from utils.enums import Likelihood, ParamTag
from utils.plot import plot_BSM_angles_limit

rc('text', usetex=False)
rc('font', **{'family':'serif', 'serif':['Computer Modern'], 'size':18})


RUN = False


z = 0.
scenarios = [
    [np.sin(np.pi/2.)**2, z, z, z],
    [z, np.cos(np.pi/2.)**4, z, z],
    [z, z, np.sin(np.pi/2.)**2, z],
]
xticks = [r'$\mathcal{O}_{12}$', r'$\mathcal{O}_{13}$', r'$\mathcal{O}_{23}$']


def main():
    args = fr.parse_args()
    fr.process_args(args)
    misc_utils.print_args(args)

    bins = 10

    asimov_paramset, mcmc_paramset = fr.get_paramsets(args)

    sc_range = mcmc_paramset.from_tag(ParamTag.SCALE)[0].ranges
    scan_scales = np.linspace(sc_range[0], sc_range[1], bins)
    print 'scan_scales', scan_scales

    p = mcmc_paramset.from_tag([ParamTag.SCALE, ParamTag.MMANGLES], invert=True)
    n_params = len(p)
    prior_ranges = p.seeds

    outfile = './sens'
    if RUN:
        if args.likelihood is Likelihood.GOLEMFIT:
            fitter = gf_utils.setup_fitter(args, asimov_paramset)
        elif args.likelihood is Likelihood.GAUSSIAN:
            fitter = None

        def CubePrior(cube, ndim, nparams):
            # default are uniform priors
            return ;

        data = np.zeros((len(scenarios), bins, 2))
        mm_angles = mcmc_paramset.from_tag(ParamTag.MMANGLES)
        sc_angles = mcmc_paramset.from_tag(ParamTag.SCALE)[0]
        for idx, scen in enumerate(scenarios):
            scales = []
            evidences = []
            for yidx, an in enumerate(mm_angles):
                an.value = scen[yidx]
            for sc in scan_scales:
                sc_angles.value = sc
                def lnProb(cube, ndim, nparams):
                    for i in range(ndim):
                        prange = prior_ranges[i][1] - prior_ranges[i][0]
                        p[i].value = prange*cube[i] + prior_ranges[i][0]
                    for name in p.names:
                        mcmc_paramset[name].value = p[name].value
                    theta = mcmc_paramset.values
                    # print 'theta', theta
                    # print 'mcmc_paramset', mcmc_paramset
                    return llh_utils.triangle_llh(
                        theta=theta,
                        args=args,
                        asimov_paramset=asimov_paramset,
                        mcmc_paramset=mcmc_paramset,
                        fitter=fitter
                    )
                # TODO(shivesh)
                prefix = 'mnrun_{0:.0E}'.format(np.power(10, sc)) + '_' + misc_utils.gen_outfile_name(args)[2:]
                print 'begin running evidence calculation for {0}'.format(prefix)
                result = pymultinest.run(
                    LogLikelihood=lnProb,
                    Prior=CubePrior,
                    n_dims=n_params,
                    importance_nested_sampling=True,
                    n_live_points=args.bayes_live_points,
                    evidence_tolerance=args.bayes_tolerance,
                    outputfiles_basename=prefix,
                    resume=False,
                    verbose=True
                )

                analyzer = pymultinest.Analyzer(outputfiles_basename=prefix, n_params=n_params)
                a_lnZ = analyzer.get_stats()['global evidence']
                print 'Evidence = {0}'.format(a_lnZ)
                scales.append(sc)
                evidences.append(a_lnZ)

            for i, d in enumerate(evidences):
                data[idx][i][0] = scales[i]
                data[idx][i][1] = d

        np.save(outfile + '.npy', data)

    plot_BSM_angles_limit(
        infile=outfile+'.npy',
        outfile=outfile,
        xticks=xticks,
        outformat=['png'],
        args=args,
        bayesian=True
    )


main.__doc__ = __doc__


if __name__ == '__main__':
    main()

