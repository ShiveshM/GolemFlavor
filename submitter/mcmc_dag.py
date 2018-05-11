#! /usr/bin/env python

import os
import numpy as np

full_scan_mfr = [
    # (1, 1, 1), (1, 0, 0)
]

fix_sfr_mfr = [
    (1, 1, 1, 1, 2, 0),
    (1, 1, 1, 1, 0, 0),
    (1, 1, 1, 0, 1, 0),
    # (1, 1, 1, 0, 0, 1),
    # (1, 1, 0, 1, 2, 0),
    # (1, 1, 0, 1, 0, 0),
    # (1, 1, 0, 0, 1, 0),
    # (1, 0, 0, 1, 0, 0),
    # (0, 1, 0, 0, 1, 0),
    # (1, 2, 0, 1, 2, 0),
    # (1, 2, 0, 0, 1, 0),
]

GLOBAL_PARAMS = {}

# MCMC
GLOBAL_PARAMS.update(dict(
    run_mcmc = 'True',
    burnin   = 250,
    nsteps   = 1000,
    nwalkers = 60,
    seed     = 25,
    mcmc_seed_type = 'uniform'
))

# FR
dimension         = [3]
# dimension         = [4, 5, 7, 8]
GLOBAL_PARAMS.update(dict(
    threads           = 1,
    binning           = '6e4 1e7 5',
    no_bsm            = 'False',
    scale_region      = "1E10",
    energy_dependance = 'spectral',
    spectral_index    = -2,
    fix_mixing        = 'False',
    fix_mixing_almost = 'False',
    fold_index        = 'True'
))

# Likelihood
GLOBAL_PARAMS.update(dict(
    likelihood  = 'golemfit',
    sigma_ratio = '0.01'
))

# GolemFit
GLOBAL_PARAMS.update(dict(
    ast  = 'p2_0',
    data = 'real'
))

# Plot
GLOBAL_PARAMS.update(dict(
    plot_angles       = 'True',
    plot_elements     = 'False',
))

outfile = 'dagman_FR_MCMC_{0}.submit'.format(GLOBAL_PARAMS['likelihood'])
golemfitsourcepath = os.environ['GOLEMSOURCEPATH'] + '/GolemFit'
condor_script = golemfitsourcepath + '/scripts/flavour_ratio/submitter/mcmc_submit.sub'

with open(outfile, 'w') as f:
    job_number = 1
    for dim in dimension:
        print 'dimension', dim
        outchain_head = '/data/user/smandalia/flavour_ratio/data/{0}/DIM{1}/'.format(
            GLOBAL_PARAMS['likelihood'], dim
        )
        for frs in fix_sfr_mfr:
            print 'frs', frs
            outchains = outchain_head + '/fix_ifr/'
            if GLOBAL_PARAMS['likelihood'].lower() == 'gaussian':
                outchains += '{0}/'.format(str(GLOBAL_PARAMS['sigma_ratio']).replace('.', '_'))
            outchains += 'mcmc_chain'
            f.write('JOB\tjob{0}\t{1}\n'.format(job_number, condor_script))
            f.write('VARS\tjob{0}\tdimension="{1}"\n'.format(job_number, dim))
            f.write('VARS\tjob{0}\tmr0="{1}"\n'.format(job_number, frs[0]))
            f.write('VARS\tjob{0}\tmr1="{1}"\n'.format(job_number, frs[1]))
            f.write('VARS\tjob{0}\tmr2="{1}"\n'.format(job_number, frs[2]))
            f.write('VARS\tjob{0}\tfix_source_ratio="{1}"\n'.format(job_number, True))
            f.write('VARS\tjob{0}\tsr0="{1}"\n'.format(job_number, frs[3]))
            f.write('VARS\tjob{0}\tsr1="{1}"\n'.format(job_number, frs[4]))
            f.write('VARS\tjob{0}\tsr2="{1}"\n'.format(job_number, frs[5]))
            for key in GLOBAL_PARAMS.iterkeys():
                f.write('VARS\tjob{0}\t{1}="{2}"\n'.format(job_number, key, GLOBAL_PARAMS[key]))
            f.write('VARS\tjob{0}\toutfile="{1}"\n'.format(job_number, outchains))
            job_number += 1

        for frs in full_scan_mfr:
            print 'frs', frs
            outchains = outchain_head + '/full/'
            if GLOBAL_PARAMS['likelihood'].lower() == 'gaussian':
                outchains += '{0}/'.format(str(GLOBAL_PARAMS['sigma_ratio']).replace('.', '_'))
            outchains += 'mcmc_chain'
            f.write('JOB\tjob{0}\t{1}\n'.format(job_number, condor_script))
            f.write('VARS\tjob{0}\tdimension="{1}"\n'.format(job_number, dim))
            f.write('VARS\tjob{0}\tmr0="{1}"\n'.format(job_number, frs[0]))
            f.write('VARS\tjob{0}\tmr1="{1}"\n'.format(job_number, frs[1]))
            f.write('VARS\tjob{0}\tmr2="{1}"\n'.format(job_number, frs[2]))
            f.write('VARS\tjob{0}\tfix_source_ratio="{1}"\n'.format(job_number, False))
            f.write('VARS\tjob{0}\tsr0="{1}"\n'.format(job_number, 0))
            f.write('VARS\tjob{0}\tsr1="{1}"\n'.format(job_number, 0))
            f.write('VARS\tjob{0}\tsr2="{1}"\n'.format(job_number, 0))
            for key in GLOBAL_PARAMS.iterkeys():
                f.write('VARS\tjob{0}\t{1}="{2}"\n'.format(job_number, key, GLOBAL_PARAMS[key]))
            f.write('VARS\tjob{0}\toutfile="{1}"\n'.format(job_number, outchains))
            job_number += 1

    print 'dag file = {0}'.format(outfile)
