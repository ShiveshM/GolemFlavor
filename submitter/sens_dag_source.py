#! /usr/bin/env python

import os
import numpy as np

full_scan_mfr = [
    # (1, 1, 1), (1, 2, 0)
]

bins = 5
binning = np.linspace(0, 1, bins)
grid = np.dstack(np.meshgrid(binning, binning)).reshape(bins*bins, 2)
sources = []
for x in grid:
    if x[0]+x[1] > 1:
        continue
    sources.append([x[0], x[1], 1-x[0]-x[1]])

# fix_sfr_mfr = [
#     (1, 1, 1, 1, 2, 0),
#     (1, 1, 1, 1, 0, 0),
#     (1, 1, 1, 0, 1, 0),
#     # (1, 1, 1, 0, 0, 1),
#     # (1, 1, 0, 1, 2, 0),
#     # (1, 1, 0, 1, 0, 0),
#     # (1, 1, 0, 0, 1, 0),
#     # (1, 0, 0, 1, 0, 0),
#     # (0, 1, 0, 0, 1, 0),
#     # (1, 2, 0, 1, 2, 0),
#     # (1, 2, 0, 0, 1, 0),
# ]
fix_sfr_mfr = []
for s in sources:
    fix_sfr_mfr.append((1, 1, 1, s[0], s[1], s[2]))
print 'fix_sfr_mfr', fix_sfr_mfr
print 'len(fix_sfr_mfr)', len(fix_sfr_mfr)

GLOBAL_PARAMS = {}

# Bayes Factor
sens_eval_bin = 'true' # set to 'all' to run normally
GLOBAL_PARAMS.update(dict(
    sens_run      = 'True',
    run_method    = 'fixed_angle', # full, fixed_angle, corr_angle
    stat_method   = 'bayesian',
    sens_bins     = 10,
    seed          = None
))

# MultiNest
GLOBAL_PARAMS.update(dict(
    # mn_live_points = 1000,
    # mn_live_points = 600,
    mn_live_points = 100,
    # mn_tolerance   = 0.1,
    mn_tolerance   = 0.3,
    mn_output      = './mnrun'
))

# FR
dimension         = [6]
# dimension         = [3, 6]
# dimension         = [3, 4, 5, 6, 7, 8]
GLOBAL_PARAMS.update(dict(
    threads            = 1,
    binning            = '6e4 1e7 20',
    no_bsm             = 'False',
    scale_region       = "1E10",
    energy_dependance  = 'spectral',
    spectral_index     = -2,
    fix_mixing         = 'None',
    fix_mixing_almost  = 'False',
    fold_index         = 'True',
    save_measured_fr   = 'False',
    output_measured_fr = './frs/'
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
    plot_statistic = 'True'
))

outfile = 'dagman_FR_SENS_{0}_{1}_{2}_{3}'.format(
    GLOBAL_PARAMS['stat_method'], GLOBAL_PARAMS['run_method'],
    GLOBAL_PARAMS['likelihood'], GLOBAL_PARAMS['data']
)
# outfile += '_seed2'
# outfile += '_tol03'
# outfile += '_NULL'
# outfile += '_prior'
# outfile += '_strictprior'
# outfile += '_noprior'
outfile += '_sourcescan'
outfile += '.submit'
golemfitsourcepath = os.environ['GOLEMSOURCEPATH'] + '/GolemFit'
condor_script = golemfitsourcepath + '/scripts/flavour_ratio/submitter/sens_submit.sub'

if sens_eval_bin.lower() != 'all':
    if GLOBAL_PARAMS['run_method'].lower() == 'corr_angle':
        raise NotImplementedError
        sens_runs = GLOBAL_PARAMS['sens_bins']**2
    else:
        sens_runs = GLOBAL_PARAMS['sens_bins'] + 1
else: sens_runs = 1

with open(outfile, 'w') as f:
    job_number = 1
    for dim in dimension:
        print 'dimension', dim
        outchain_head = '/data/user/smandalia/flavour_ratio/data/{0}/DIM{1}'.format(
            GLOBAL_PARAMS['likelihood'], dim
        )
        for frs in fix_sfr_mfr:
            print 'frs', frs
            output = outchain_head + '/fix_ifr/'
            if GLOBAL_PARAMS['likelihood'].lower() == 'gaussian':
                output += '{0}/'.format(str(GLOBAL_PARAMS['sigma_ratio']).replace('.', '_'))
            # output += 'seed2/'
            # output += 'mn_noverlap/'
            # output += 'tol_03/'
            # output += 'prior/'
            # output += 'strictprior/'
            # output += 'noprior/'
            output += 'sourcescan/'
            for r in xrange(sens_runs):
                print 'run', r
                f.write('JOB\tjob{0}\t{1}\n'.format(job_number, condor_script))
                f.write('VARS\tjob{0}\tdimension="{1}"\n'.format(job_number, dim))
                f.write('VARS\tjob{0}\tmr0="{1}"\n'.format(job_number, frs[0]))
                f.write('VARS\tjob{0}\tmr1="{1}"\n'.format(job_number, frs[1]))
                f.write('VARS\tjob{0}\tmr2="{1}"\n'.format(job_number, frs[2]))
                f.write('VARS\tjob{0}\tfix_source_ratio="{1}"\n'.format(job_number, True))
                f.write('VARS\tjob{0}\tsr0="{1}"\n'.format(job_number, frs[3]))
                f.write('VARS\tjob{0}\tsr1="{1}"\n'.format(job_number, frs[4]))
                f.write('VARS\tjob{0}\tsr2="{1}"\n'.format(job_number, frs[5]))
                if sens_eval_bin.lower() != 'all':
                    f.write('VARS\tjob{0}\tsens_eval_bin="{1}"\n'.format(job_number, r))
                else:
                    f.write('VARS\tjob{0}\tsens_eval_bin="{1}"\n'.format(job_number, 'all'))
                for key in GLOBAL_PARAMS.iterkeys():
                    f.write('VARS\tjob{0}\t{1}="{2}"\n'.format(job_number, key, GLOBAL_PARAMS[key]))
                f.write('VARS\tjob{0}\toutfile="{1}"\n'.format(job_number, output))
                job_number += 1
                # break

        # for frs in full_scan_mfr:
        #     print 'frs', frs
        #     output = outchain_head + '/full/'
        #     if GLOBAL_PARAMS['likelihood'].lower() == 'gaussian':
        #         output += '{0}/'.format(str(GLOBAL_PARAMS['sigma_ratio']).replace('.', '_'))
        #     for r in xrange(sens_runs):
        #         print 'run', r
        #         f.write('JOB\tjob{0}\t{1}\n'.format(job_number, condor_script))
        #         f.write('VARS\tjob{0}\tdimension="{1}"\n'.format(job_number, dim))
        #         f.write('VARS\tjob{0}\tmr0="{1}"\n'.format(job_number, frs[0]))
        #         f.write('VARS\tjob{0}\tmr1="{1}"\n'.format(job_number, frs[1]))
        #         f.write('VARS\tjob{0}\tmr2="{1}"\n'.format(job_number, frs[2]))
        #         f.write('VARS\tjob{0}\tfix_source_ratio="{1}"\n'.format(job_number, False))
        #         f.write('VARS\tjob{0}\tsr0="{1}"\n'.format(job_number, 0))
        #         f.write('VARS\tjob{0}\tsr1="{1}"\n'.format(job_number, 0))
        #         f.write('VARS\tjob{0}\tsr2="{1}"\n'.format(job_number, 0))
        #         if sens_eval_bin.lower() != 'all':
        #             f.write('VARS\tjob{0}\tsens_eval_bin="{1}"\n'.format(job_number, r))
        #         else:
        #             f.write('VARS\tjob{0}\tsens_eval_bin="{1}"\n'.format(job_number, 'all'))
        #         for key in GLOBAL_PARAMS.iterkeys():
        #             f.write('VARS\tjob{0}\t{1}="{2}"\n'.format(job_number, key, GLOBAL_PARAMS[key]))
        #         f.write('VARS\tjob{0}\toutfile="{1}"\n'.format(job_number, output))
        #         job_number += 1

    print 'dag file = {0}'.format(outfile)