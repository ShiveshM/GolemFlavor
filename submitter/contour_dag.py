#! /usr/bin/env python

import os
import numpy as np

gfsource      = os.environ['GOLEMSOURCEPATH'] + '/GolemFit'
condor_script = gfsource + '/scripts/flavour_ratio/submitter/contour_submit.sub'

injected_ratios = [
    (1, 1, 1),
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1)
]

GLOBAL_PARAMS = {}

GLOBAL_PARAMS.update(dict(
    threads            = 1,
    save_measured_fr   = 'False',
    output_measured_fr = './frs/',
    seed               = None
))

# MultiNest
GLOBAL_PARAMS.update(dict(
    mn_live_points = 5000,
    mn_tolerance   = 0.3,
))

# Likelihood
GLOBAL_PARAMS.update(dict(
    likelihood  = 'golemfit',
))

# GolemFit
GLOBAL_PARAMS.update(dict(
    ast  = 'p2_0',
    # data = 'realisation'
    # data = 'asimov'
    data = 'real'
))

# Plot
GLOBAL_PARAMS.update(dict(
    plot_chains   = 'False',
    plot_triangle = 'False'
))

outfile = 'dagman_FR_CONTOUR_{0}'.format(GLOBAL_PARAMS['data'])
outfile += '.submit'
output  = '/data/user/smandalia/flavour_ratio/data/contour/{0}/{1}/'.format(
    GLOBAL_PARAMS['likelihood'], GLOBAL_PARAMS['data']
)
# output += 'nosyst/'
# output += 'noprompt/'
# output += 'strictpriors/'

with open(outfile, 'w') as f:
    job_number = 1
    for inj in injected_ratios:
        print 'inj', inj
        f.write('JOB\tjob{0}\t{1}\n'.format(job_number, condor_script))
        f.write('VARS\tjob{0}\tir0="{1}"\n'.format(job_number, inj[0]))
        f.write('VARS\tjob{0}\tir1="{1}"\n'.format(job_number, inj[1]))
        f.write('VARS\tjob{0}\tir2="{1}"\n'.format(job_number, inj[2]))
        for key in GLOBAL_PARAMS.iterkeys():
            f.write('VARS\tjob{0}\t{1}="{2}"\n'.format(job_number, key, GLOBAL_PARAMS[key]))
        f.write('VARS\tjob{0}\toutfile="{1}"\n'.format(job_number, output))
        job_number += 1
        if GLOBAL_PARAMS['data'] == 'real': break

print 'dag file = {0}'.format(outfile)
