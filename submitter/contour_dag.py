#! /usr/bin/env python

from __future__ import absolute_import, division, print_function

import os
import numpy as np

injected_ratios = [
    (1, 1, 1),
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1)
]

datadir = '/data/user/smandalia/flavour_ratio/data/contour'

prefix = ''
# prefix = '_noprior'

golemfitsourcepath = os.environ['GOLEMSOURCEPATH'] + '/GolemFit'
condor_script = golemfitsourcepath + '/scripts/flavour_ratio/submitter/contour_submit.sub'

GLOBAL_PARAMS = {}

GLOBAL_PARAMS.update(dict(
    threads     = 12,
    seed        = 26
))

# Emcee
GLOBAL_PARAMS.update(dict(
    run_mcmc = 'True',
    burnin   = 1500,
    nsteps   = 15000,
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

dagfile = 'dagman_CONTOUR_{0}'.format(GLOBAL_PARAMS['data'])
dagfile += prefix + '.submit'

with open(dagfile, 'w') as f:
    job_number = 1
    for inj in injected_ratios:
        print('inj', inj)
        f.write('JOB\tjob{0}\t{1}\n'.format(job_number, condor_script))
        f.write('VARS\tjob{0}\tir0="{1}"\n'.format(job_number, inj[0]))
        f.write('VARS\tjob{0}\tir1="{1}"\n'.format(job_number, inj[1]))
        f.write('VARS\tjob{0}\tir2="{1}"\n'.format(job_number, inj[2]))
        for key in GLOBAL_PARAMS.iterkeys():
            f.write('VARS\tjob{0}\t{1}="{2}"\n'.format(job_number, key, GLOBAL_PARAMS[key]))
        f.write('VARS\tjob{0}\tdatadir="{1}"\n'.format(job_number, datadir))
        job_number += 1
        if GLOBAL_PARAMS['data'] == 'real': break

print('total jobs = {0}'.format(job_number - 1))
print('dag file = {0}'.format(dagfile))
