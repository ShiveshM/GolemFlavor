#! /usr/bin/env python

import os
import numpy as np

a_fr = (1, 2, 0)
b_fr = (1, 0, 0)
c_fr = (0, 1, 0)
d_fr = (0, 0, 1)
e_fr = (1, 1, 1)
f_fr = (2, 1, 0)
g_fr = (1, 1, 0)

full_scan_mfr = [
    (1, 1, 1), (1, 1, 0)
]

fix_sfr_mfr = [
    (1, 1, 1, 1, 0, 0),
    (1, 1, 1, 0, 1, 0),
    # (1, 1, 1, 0, 0, 1),
    (1, 1, 1, 1, 2, 0),
    # (1, 1, 0, 0, 1, 0),
    (1, 1, 0, 1, 2, 0),
    # (1, 1, 0, 1, 0, 0),
    # (1, 0, 0, 1, 0, 0),
    (0, 1, 0, 0, 1, 0),
    # (1, 2, 0, 0, 1, 0),
    # (1, 2, 0, 1, 2, 0)
]

# MCMC
run_mcmc = 'True'
burnin   = 1000
nsteps   = 4000
nwalkers = 70
seed     = 24
threads  = 12
mcmc_seed_type = 'uniform'

# FR
dimension         = [3, 6]
energy            = [1e6]
likelihood        = 'golemfit'
no_bsm            = 'False'
sigma_ratio       = ['0.01']
scale             = "1E-20 1E-30"
scale_region      = "1E10"
energy_dependance = 'spectral'
spectral_index    = -2
binning           = [4, 7, 51]

# Likelihood
likelihood = 'golemfit'

# Nuisance
astroDeltaGamma = 2.
astroNorm       = 1.
convNorm        = 1.
muonNorm        = 1.
promptNorm      = 0.

# GolemFit
ast  = 'p2_0'
data = 'real'

# Plot
plot_angles   = 'True'
plot_elements = 'False'

outfile = 'dagman_FR.submit'
golemfitsourcepath = os.environ['GOLEMSOURCEPATH'] + '/GolemFit'
condor_script = golemfitsourcepath + '/scripts/flavour_ratio/submitter/submit.sub'

with open(outfile, 'w') as f:
    job_number = 1
    for dim in dimension:
        print 'dimension', dim
        for en in energy:
            print 'energy {0:.0E}'.format(en)

            outchain_head = '/data/user/smandalia/flavour_ratio/data/{0}/DIM{1}/{2:.0E}'.format(likelihood, dim, en)

            for sig in sigma_ratio:
                print 'sigma', sig
                for frs in fix_sfr_mfr:
                    print frs
                    outchains = outchain_head + '/fix_ifr/{0}/mcmc_chain'.format(str(sig).replace('.', '_'))
                    f.write('JOB\tjob{0}\t{1}\n'.format(job_number, condor_script))
                    f.write('VARS\tjob{0}\tmr0="{1}"\n'.format(job_number, frs[0]))
                    f.write('VARS\tjob{0}\tmr1="{1}"\n'.format(job_number, frs[1]))
                    f.write('VARS\tjob{0}\tmr2="{1}"\n'.format(job_number, frs[2]))
                    f.write('VARS\tjob{0}\tsigma_ratio="{1}"\n'.format(job_number, sig))
                    f.write('VARS\tjob{0}\tfix_source_ratio="{1}"\n'.format(job_number, 'True'))
                    f.write('VARS\tjob{0}\tsr0="{1}"\n'.format(job_number, frs[3]))
                    f.write('VARS\tjob{0}\tsr1="{1}"\n'.format(job_number, frs[4]))
                    f.write('VARS\tjob{0}\tsr2="{1}"\n'.format(job_number, frs[5]))
                    f.write('VARS\tjob{0}\tfix_scale="{1}"\n'.format(job_number, 'False'))
                    f.write('VARS\tjob{0}\tscale="{1}"\n'.format(job_number, 0))
                    f.write('VARS\tjob{0}\tscale_region="{1}"\n'.format(job_number, scale_region))
                    f.write('VARS\tjob{0}\tdimension="{1}"\n'.format(job_number, dim))
                    f.write('VARS\tjob{0}\tenergy="{1}"\n'.format(job_number, en))
                    f.write('VARS\tjob{0}\tlikelihood="{1}"\n'.format(job_number, likelihood))
                    f.write('VARS\tjob{0}\tburnin="{1}"\n'.format(job_number, burnin))
                    f.write('VARS\tjob{0}\tnwalkers="{1}"\n'.format(job_number, nwalkers))
                    f.write('VARS\tjob{0}\tnsteps="{1}"\n'.format(job_number, nsteps))
                    f.write('VARS\tjob{0}\toutfile="{1}"\n'.format(job_number, outchains))
                    f.write('VARS\tjob{0}\tfix_mixing="{1}"\n'.format(job_number, 'False'))
                    f.write('VARS\tjob{0}\tno_bsm="{1}"\n'.format(job_number, no_bsm))
                    f.write('VARS\tjob{0}\trun_mcmc="{1}"\n'.format(job_number, run_mcmc))
                    f.write('VARS\tjob{0}\tastroDeltaGamma="{1}"\n'.format(job_number, astroDeltaGamma))
                    f.write('VARS\tjob{0}\tastroNorm="{1}"\n'.format(job_number, astroNorm))
                    f.write('VARS\tjob{0}\tconvNorm="{1}"\n'.format(job_number, convNorm))
                    f.write('VARS\tjob{0}\tmuonNorm="{1}"\n'.format(job_number, muonNorm))
                    f.write('VARS\tjob{0}\tpromptNorm="{1}"\n'.format(job_number, promptNorm))
                    f.write('VARS\tjob{0}\tdata="{1}"\n'.format(job_number, data))
                    f.write('VARS\tjob{0}\tast="{1}"\n'.format(job_number, ast))
                    f.write('VARS\tjob{0}\tplot_angles="{1}"\n'.format(job_number, plot_angles))
                    f.write('VARS\tjob{0}\tplot_elements="{1}"\n'.format(job_number, plot_elements))
                    f.write('VARS\tjob{0}\tseed="{1}"\n'.format(job_number, seed))
                    f.write('VARS\tjob{0}\tthreads="{1}"\n'.format(job_number, threads))
                    f.write('VARS\tjob{0}\tlikelihood="{1}"\n'.format(job_number, likelihood))
                    f.write('VARS\tjob{0}\tmcmc_seed_type="{1}"\n'.format(job_number, mcmc_seed_type))
                    f.write('VARS\tjob{0}\tenergy_dependance="{1}"\n'.format(job_number, energy_dependance))
                    f.write('VARS\tjob{0}\tspectral_index="{1}"\n'.format(job_number, spectral_index))
                    f.write('VARS\tjob{0}\tbinning_0="{1}"\n'.format(job_number, binning[0]))
                    f.write('VARS\tjob{0}\tbinning_1="{1}"\n'.format(job_number, binning[1]))
                    f.write('VARS\tjob{0}\tbinning_2="{1}"\n'.format(job_number, binning[2]))
                    job_number += 1

                for frs in full_scan_mfr:
                    print frs
                    outchains = outchain_head + '/full_scan/{0}/mcmc_chain'.format(str(sig).replace('.', '_'))
                    f.write('JOB\tjob{0}\t{1}\n'.format(job_number, condor_script))
                    f.write('VARS\tjob{0}\tmr0="{1}"\n'.format(job_number, frs[0]))
                    f.write('VARS\tjob{0}\tmr1="{1}"\n'.format(job_number, frs[1]))
                    f.write('VARS\tjob{0}\tmr2="{1}"\n'.format(job_number, frs[2]))
                    f.write('VARS\tjob{0}\tsigma_ratio="{1}"\n'.format(job_number, sig))
                    f.write('VARS\tjob{0}\tfix_source_ratio="{1}"\n'.format(job_number, 'False'))
                    f.write('VARS\tjob{0}\tsr0="{1}"\n'.format(job_number, 0))
                    f.write('VARS\tjob{0}\tsr1="{1}"\n'.format(job_number, 0))
                    f.write('VARS\tjob{0}\tsr2="{1}"\n'.format(job_number, 0))
                    f.write('VARS\tjob{0}\tfix_scale="{1}"\n'.format(job_number, 'False'))
                    f.write('VARS\tjob{0}\tscale="{1}"\n'.format(job_number, 0))
                    f.write('VARS\tjob{0}\tscale_region="{1}"\n'.format(job_number, scale_region))
                    f.write('VARS\tjob{0}\tdimension="{1}"\n'.format(job_number, dim))
                    f.write('VARS\tjob{0}\tenergy="{1}"\n'.format(job_number, en))
                    f.write('VARS\tjob{0}\tlikelihood="{1}"\n'.format(job_number, likelihood))
                    f.write('VARS\tjob{0}\tburnin="{1}"\n'.format(job_number, burnin))
                    f.write('VARS\tjob{0}\tnwalkers="{1}"\n'.format(job_number, nwalkers))
                    f.write('VARS\tjob{0}\tnsteps="{1}"\n'.format(job_number, nsteps))
                    f.write('VARS\tjob{0}\toutfile="{1}"\n'.format(job_number, outchains))
                    f.write('VARS\tjob{0}\tfix_mixing="{1}"\n'.format(job_number, 'False'))
                    f.write('VARS\tjob{0}\tno_bsm="{1}"\n'.format(job_number, no_bsm))
                    f.write('VARS\tjob{0}\trun_mcmc="{1}"\n'.format(job_number, run_mcmc))
                    f.write('VARS\tjob{0}\tastroDeltaGamma="{1}"\n'.format(job_number, astroDeltaGamma))
                    f.write('VARS\tjob{0}\tastroNorm="{1}"\n'.format(job_number, astroNorm))
                    f.write('VARS\tjob{0}\tconvNorm="{1}"\n'.format(job_number, convNorm))
                    f.write('VARS\tjob{0}\tmuonNorm="{1}"\n'.format(job_number, muonNorm))
                    f.write('VARS\tjob{0}\tpromptNorm="{1}"\n'.format(job_number, promptNorm))
                    f.write('VARS\tjob{0}\tdata="{1}"\n'.format(job_number, data))
                    f.write('VARS\tjob{0}\tast="{1}"\n'.format(job_number, ast))
                    f.write('VARS\tjob{0}\tplot_angles="{1}"\n'.format(job_number, plot_angles))
                    f.write('VARS\tjob{0}\tplot_elements="{1}"\n'.format(job_number, plot_elements))
                    f.write('VARS\tjob{0}\tseed="{1}"\n'.format(job_number, seed))
                    f.write('VARS\tjob{0}\tthreads="{1}"\n'.format(job_number, threads))
                    f.write('VARS\tjob{0}\tlikelihood="{1}"\n'.format(job_number, likelihood))
                    f.write('VARS\tjob{0}\tmcmc_seed_type="{1}"\n'.format(job_number, mcmc_seed_type))
                    f.write('VARS\tjob{0}\tenergy_dependance="{1}"\n'.format(job_number, energy_dependance))
                    f.write('VARS\tjob{0}\tspectral_index="{1}"\n'.format(job_number, spectral_index))
                    f.write('VARS\tjob{0}\tbinning_0="{1}"\n'.format(job_number, binning[0]))
                    f.write('VARS\tjob{0}\tbinning_1="{1}"\n'.format(job_number, binning[1]))
                    f.write('VARS\tjob{0}\tbinning_2="{1}"\n'.format(job_number, binning[2]))
                    job_number += 1
