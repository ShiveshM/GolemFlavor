#! /usr/bin/env python

import os
import numpy as np

scenarios = [
    ((0, 1, 0), 'OET'),
    ((1, 0, 0), 'OUT')
]

dims = [
    3, 4, 5, 7, 8
]

datadir = '/data/user/smandalia/flavour_ratio/data/fr'

prefix = ''
# prefix = '_noprior'

golemfitsourcepath = os.environ['GOLEMSOURCEPATH'] + '/GolemFit'
condor_script = golemfitsourcepath + '/scripts/flavour_ratio/submitter/fr_submit.sub'

GLOBAL_PARAMS = {}

# General
GLOBAL_PARAMS.update(dict(
    threads     = 12,
    seed        = 26
))

# FR
GLOBAL_PARAMS.update(dict(
    binning = '6e4 1e7 20',
    no_bsm  = 'False'
))

# Emcee
GLOBAL_PARAMS.update(dict(
    run_mcmc = 'True',
    burnin   = 1000,
    nsteps   = 10000,
    nwalkers = 60,
    mcmc_seed_type = 'uniform'
))

# GolemFit
GLOBAL_PARAMS.update(dict(
    ast  = 'p2_0',
    data = 'real'
))

# Plot
GLOBAL_PARAMS.update(dict(
    plot_angles   = 'False',
    plot_elements = 'False',
))

dagfile = 'dagman_FR_{0}'.format(GLOBAL_PARAMS['data'])
dagfile += prefix + '.submit'

with open(dagfile, 'w') as f:
    job_number = 1
    for dim in dims:
        print 'dims', dim
        of_d = datadir + '/DIM{0}/{1}'.format(dim, prefix)
        for src, tex in scenarios:
            print 'scenario: src =', src, 'tex =', tex
            f.write('JOB\tjob{0}\t{1}\n'.format(job_number, condor_script))
            f.write('VARS\tjob{0}\tdimension="{1}"\n'.format(job_number, dim))
            f.write('VARS\tjob{0}\ttexture="{1}"\n'.format(job_number, tex))
            f.write('VARS\tjob{0}\tsr0="{1}"\n'.format(job_number, src[0]))
            f.write('VARS\tjob{0}\tsr1="{1}"\n'.format(job_number, src[1]))
            f.write('VARS\tjob{0}\tsr2="{1}"\n'.format(job_number, src[2]))
            for key in GLOBAL_PARAMS.iterkeys():
                f.write('VARS\tjob{0}\t{1}="{2}"\n'.format(job_number, key, GLOBAL_PARAMS[key]))
            f.write('VARS\tjob{0}\tdatadir="{1}"\n'.format(job_number, datadir))
            job_number += 1

print 'total jobs = {0}'.format(job_number - 1)
print 'dag file = {0}'.format(dagfile)

