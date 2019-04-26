#! /usr/bin/env python

import os
import numpy as np

source_ratios = [
    (1, 0, 0),
    (0, 1, 0)
]

textures = [
    'OET', 'OUT'
]

datadir = '/data/user/smandalia/flavour_ratio/data/mc_texture'

prefix = ''
# prefix = '_noprior'

golemfitsourcepath = os.environ['GOLEMSOURCEPATH'] + '/GolemFit'
condor_script = golemfitsourcepath + '/scripts/flavour_ratio/submitter/mc_texture_submit.sub'

GLOBAL_PARAMS = {}

GLOBAL_PARAMS.update(dict(
    threads     = 1,
    seed        = 26
))

# Emcee
GLOBAL_PARAMS.update(dict(
    run_mcmc = 'True',
    burnin   = 200,
    nsteps   = 1000,
    nwalkers = 60,
    mcmc_seed_type = 'uniform'
))

# FR
GLOBAL_PARAMS.update(dict(
    binning    = '6e4 1e7 20',
    dimension  = 6,
    no_bsm     = 'False'
))

# Plot
GLOBAL_PARAMS.update(dict(
    plot_angles   = 'False',
    plot_elements = 'False',
))

dagfile = 'dagman_mc_texture'
dagfile += prefix + '.submit'

with open(dagfile, 'w') as f:
    job_number = 1
    for src in source_ratios:
        print 'src', src
        for tex in textures:
            print 'texture', tex
            f.write('JOB\tjob{0}\t{1}\n'.format(job_number, condor_script))
            f.write('VARS\tjob{0}\tsr0="{1}"\n'.format(job_number, src[0]))
            f.write('VARS\tjob{0}\tsr1="{1}"\n'.format(job_number, src[1]))
            f.write('VARS\tjob{0}\tsr2="{1}"\n'.format(job_number, src[2]))
            f.write('VARS\tjob{0}\ttexture="{1}"\n'.format(job_number, tex))
            for key in GLOBAL_PARAMS.iterkeys():
                f.write('VARS\tjob{0}\t{1}="{2}"\n'.format(job_number, key, GLOBAL_PARAMS[key]))
            f.write('VARS\tjob{0}\tdatadir="{1}"\n'.format(job_number, datadir))
            job_number += 1

print 'total jobs = {0}'.format(job_number - 1)
print 'dag file = {0}'.format(dagfile)

