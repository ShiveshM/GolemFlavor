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
    # (1, 1, 1), (1, 1, 0)
]

fix_sfr_mfr = [
    (1, 1, 1, 1, 2, 0),
    (1, 1, 1, 1, 0, 0),
    (1, 1, 1, 0, 1, 0),
    (1, 1, 1, 0, 0, 1),
    # (1, 1, 0, 1, 2, 0),
    # (1, 1, 0, 1, 0, 0),
    # (1, 1, 0, 0, 1, 0),
    # (1, 0, 0, 1, 0, 0),
    # (0, 1, 0, 0, 1, 0),
    # (1, 2, 0, 1, 2, 0),
    # (1, 2, 0, 0, 1, 0),
]

# MCMC
run_mcmc = 'True'
burnin   = 2500
nsteps   = 10000
nwalkers = 60
seed     = 'None'
threads  = 1
mcmc_seed_type = 'uniform'

# FR
dimension         = [3, 6]
energy            = [1e6]
no_bsm            = 'False'
sigma_ratio       = ['0.01']
scale             = "1E-20 1E-30"
scale_region      = "1E10"
energy_dependance = 'spectral'
spectral_index    = -2
binning           = [1e4, 1e7, 5]
fix_mixing        = 'False'
fix_mixing_almost = 'False'

# Likelihood
likelihood = 'gaussian'

# Nuisance
convNorm        = 1.
promptNorm      = 0.
muonNorm        = 1.
astroNorm       = 6.9
astroDeltaGamma = 2.5

# GolemFit
ast  = 'p2_0'
data = 'real'

# Bayes Factor
run_bayes_factor       = 'False'
run_angles_limit       = 'False'
run_angles_correlation = 'False'
bayes_bins             = 100
bayes_live_points      = 3000
bayes_tolerance        = 0.01
bayes_eval_bin         = 'all' # set to 'all' to run normally

# Plot
plot_angles       = 'True'
plot_elements     = 'False'
plot_bayes        = 'False'
plot_angles_limit = 'False'

outfile = 'dagman_FR.submit'
golemfitsourcepath = os.environ['GOLEMSOURCEPATH'] + '/GolemFit'
condor_script = golemfitsourcepath + '/scripts/flavour_ratio/submitter/submit.sub'

if bayes_eval_bin != 'all':
    if run_angles_correlation == 'True':
        b_runs = bayes_bins**2
    else:
        b_runs = bayes_bins
else: b_runs = 1

with open(outfile, 'w') as f:
    job_number = 1
    for dim in dimension:
        print 'dimension', dim
        for en in energy:
            print 'energy {0:.0E}'.format(en)

            if energy_dependance == 'mono':
                outchain_head = '/data/user/smandalia/flavour_ratio/data/{0}/DIM{1}/{2:.0E}'.format(likelihood, dim, en)
            elif energy_dependance == 'spectral':
                outchain_head = '/data/user/smandalia/flavour_ratio/data/{0}/DIM{1}/SI_{2}'.format(likelihood, dim, spectral_index)

            bayes_output = 'None'
            angles_lim_output = 'None'
            angles_corr_output = 'None'
            for sig in sigma_ratio:
                print 'sigma', sig
                for frs in fix_sfr_mfr:
                    print frs
                    outchains = outchain_head + '/fix_ifr/{0}/'.format(str(sig).replace('.', '_'))
                    if run_bayes_factor == 'True':
                        bayes_output = outchains + '/bayes_factor/'
                    if run_angles_limit == 'True':
                        angles_lim_output = outchains + '/angles_limit/'
                    if run_angles_correlation == 'True':
                        angles_corr_output = outchains + '/angles_corr/'
                    outchains += 'mcmc_chain'
                    for r in range(b_runs):
                        print 'run', r
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
                        f.write('VARS\tjob{0}\tfix_mixing="{1}"\n'.format(job_number, fix_mixing))
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
                        f.write('VARS\tjob{0}\tmcmc_seed_type="{1}"\n'.format(job_number, mcmc_seed_type))
                        f.write('VARS\tjob{0}\tenergy_dependance="{1}"\n'.format(job_number, energy_dependance))
                        f.write('VARS\tjob{0}\tspectral_index="{1}"\n'.format(job_number, spectral_index))
                        f.write('VARS\tjob{0}\tbinning_0="{1}"\n'.format(job_number, binning[0]))
                        f.write('VARS\tjob{0}\tbinning_1="{1}"\n'.format(job_number, binning[1]))
                        f.write('VARS\tjob{0}\tbinning_2="{1}"\n'.format(job_number, binning[2]))
                        f.write('VARS\tjob{0}\tfix_mixing_almost="{1}"\n'.format(job_number, fix_mixing_almost))
                        # f.write('VARS\tjob{0}\trun_bayes_factor="{1}"\n'.format(job_number, run_bayes_factor))
                        # f.write('VARS\tjob{0}\tbayes_bins="{1}"\n'.format(job_number, bayes_bins))
                        # f.write('VARS\tjob{0}\tbayes_output="{1}"\n'.format(job_number, bayes_output))
                        # f.write('VARS\tjob{0}\tbayes_live_points="{1}"\n'.format(job_number, bayes_live_points))
                        # f.write('VARS\tjob{0}\tbayes_tolerance="{1}"\n'.format(job_number, bayes_tolerance))
                        # f.write('VARS\tjob{0}\tplot_bayes="{1}"\n'.format(job_number, plot_bayes))
                        # f.write('VARS\tjob{0}\tbayes_eval_bin="{1}"\n'.format(job_number, r))
                        # f.write('VARS\tjob{0}\trun_angles_limit="{1}"\n'.format(job_number, run_angles_limit))
                        # f.write('VARS\tjob{0}\tangles_lim_output="{1}"\n'.format(job_number, angles_lim_output))
                        # f.write('VARS\tjob{0}\tplot_angles_limit="{1}"\n'.format(job_number, plot_angles_limit))
                        # f.write('VARS\tjob{0}\trun_angles_correlation="{1}"\n'.format(job_number, run_angles_correlation))
                        # f.write('VARS\tjob{0}\tangles_corr_output="{1}"\n'.format(job_number, angles_corr_output))
                        job_number += 1

                for frs in full_scan_mfr:
                    print frs
                    outchains = outchain_head + '/full_scan/{0}'.format(str(sig).replace('.', '_'))
                    if run_bayes_factor == 'True':
                        bayes_output = outchains + '/bayes_factor/'
                    if run_angles_limit == 'True':
                        angles_lim_output = outchains + '/angles_limit/'
                    if run_angles_correlation == 'True':
                        angles_corr_output = outchains + '/angles_corr/'
                    outchains += 'mcmc_chain'
                    for r in range(b_runs):
                        print 'run', r
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
                        f.write('VARS\tjob{0}\tfix_mixing="{1}"\n'.format(job_number, fix_mixing))
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
                        f.write('VARS\tjob{0}\tmcmc_seed_type="{1}"\n'.format(job_number, mcmc_seed_type))
                        f.write('VARS\tjob{0}\tenergy_dependance="{1}"\n'.format(job_number, energy_dependance))
                        f.write('VARS\tjob{0}\tspectral_index="{1}"\n'.format(job_number, spectral_index))
                        f.write('VARS\tjob{0}\tbinning_0="{1}"\n'.format(job_number, binning[0]))
                        f.write('VARS\tjob{0}\tbinning_1="{1}"\n'.format(job_number, binning[1]))
                        f.write('VARS\tjob{0}\tbinning_2="{1}"\n'.format(job_number, binning[2]))
                        f.write('VARS\tjob{0}\tfix_mixing_almost="{1}"\n'.format(job_number, fix_mixing_almost))
                        f.write('VARS\tjob{0}\trun_bayes_factor="{1}"\n'.format(job_number, run_bayes_factor))
                        f.write('VARS\tjob{0}\tbayes_bins="{1}"\n'.format(job_number, bayes_bins))
                        f.write('VARS\tjob{0}\tbayes_output="{1}"\n'.format(job_number, bayes_output))
                        f.write('VARS\tjob{0}\tbayes_live_points="{1}"\n'.format(job_number, bayes_live_points))
                        f.write('VARS\tjob{0}\tbayes_tolerance="{1}"\n'.format(job_number, bayes_tolerance))
                        f.write('VARS\tjob{0}\tplot_bayes="{1}"\n'.format(job_number, plot_bayes))
                        f.write('VARS\tjob{0}\tbayes_eval_bin="{1}"\n'.format(job_number, r))
                        f.write('VARS\tjob{0}\trun_angles_limit="{1}"\n'.format(job_number, run_angles_limit))
                        f.write('VARS\tjob{0}\tangles_lim_output="{1}"\n'.format(job_number, angles_lim_output))
                        f.write('VARS\tjob{0}\tplot_angles_limit="{1}"\n'.format(job_number, plot_angles_limit))
                        f.write('VARS\tjob{0}\trun_angles_correlation="{1}"\n'.format(job_number, run_angles_correlation))
                        f.write('VARS\tjob{0}\tangles_corr_output="{1}"\n'.format(job_number, angles_corr_output))
                        job_number += 1
