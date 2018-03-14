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
    # (1, 1, 1, 1, 0, 0),
    # (1, 1, 1, 0, 1, 0),
    # (1, 1, 1, 0, 0, 1),
    (1, 1, 1, 1, 2, 0),
    # (1, 1, 0, 0, 1, 0),
    # (1, 1, 0, 1, 2, 0),
    # (1, 1, 0, 1, 0, 0),
    # (1, 0, 0, 1, 0, 0),
    # (0, 1, 0, 0, 1, 0),
    # (1, 2, 0, 0, 1, 0),
    # (1, 2, 0, 1, 2, 0)
]

sigmas = ['0.01']
dimensions = [3]
energy = [1e6]
flat = False
burnin = 20
nwalkers = 200
nsteps = 100
scales = "1E-20 1E-30"
no_bsm = False

outfile = 'dagman_FR.submit'
golemfitsourcepath = os.environ['GOLEMSOURCEPATH'] + '/GolemFit'
condor_script = golemfitsourcepath + '/scripts/flavour_ratio/submitter/submit.sub'

with open(outfile, 'w') as f:
    job_number = 1
    for dim in dimensions:
        print 'dimension', dim
        for en in energy:
            print 'energy {0:.0E}'.format(en)

            outchain_head = '/data/user/smandalia/flavour_ratio/data/DIM{0}/{1:.0E}'.format(dim, en)

            for sig in sigmas:
                print 'sigma', sig
                for frs in fix_sfr_mfr:
                    print frs
                    outchains = outchain_head + '/fix_ifr/{0}/mcmc_chain'.format(str(sig).replace('.', '_'))
                    f.write('JOB\tjob{0}\t{1}\n'.format(job_number, condor_script))
                    f.write('VARS\tjob{0}\tmr0="{1}"\n'.format(job_number, frs[0]))
                    f.write('VARS\tjob{0}\tmr1="{1}"\n'.format(job_number, frs[1]))
                    f.write('VARS\tjob{0}\tmr2="{1}"\n'.format(job_number, frs[2]))
                    f.write('VARS\tjob{0}\tsigma="{1}"\n'.format(job_number, sig))
                    f.write('VARS\tjob{0}\tfix_source_ratio="{1}"\n'.format(job_number, 'True'))
                    f.write('VARS\tjob{0}\tsr0="{1}"\n'.format(job_number, frs[3]))
                    f.write('VARS\tjob{0}\tsr1="{1}"\n'.format(job_number, frs[4]))
                    f.write('VARS\tjob{0}\tsr2="{1}"\n'.format(job_number, frs[5]))
                    f.write('VARS\tjob{0}\tfix_scale="{1}"\n'.format(job_number, 'False'))
                    f.write('VARS\tjob{0}\tscale="{1}"\n'.format(job_number, 0))
                    f.write('VARS\tjob{0}\tdimension="{1}"\n'.format(job_number, dim))
                    f.write('VARS\tjob{0}\tenergy="{1}"\n'.format(job_number, en))
                    f.write('VARS\tjob{0}\tflat_llh="{1}"\n'.format(job_number, flat))
                    f.write('VARS\tjob{0}\tburnin="{1}"\n'.format(job_number, burnin))
                    f.write('VARS\tjob{0}\tnwalkers="{1}"\n'.format(job_number, nwalkers))
                    f.write('VARS\tjob{0}\tnsteps="{1}"\n'.format(job_number, nsteps))
                    f.write('VARS\tjob{0}\toutfile="{1}"\n'.format(job_number, outchains))
                    f.write('VARS\tjob{0}\tfix_mixing="{1}"\n'.format(job_number, 'False'))
                    f.write('VARS\tjob{0}\tno_bsm="{1}"\n'.format(job_number, no_bsm))
                    job_number += 1

                # for frs in full_scan_mfr:
                #     print frs
                #     outchains = outchain_head + '/full_scan/{0}/mcmc_chain'.format(str(sig).replace('.', '_'))
                #     f.write('JOB\tjob{0}\t{1}\n'.format(job_number, condor_script))
                #     f.write('VARS\tjob{0}\tmr0="{1}"\n'.format(job_number, frs[0]))
                #     f.write('VARS\tjob{0}\tmr1="{1}"\n'.format(job_number, frs[1]))
                #     f.write('VARS\tjob{0}\tmr2="{1}"\n'.format(job_number, frs[2]))
                #     f.write('VARS\tjob{0}\tsigma="{1}"\n'.format(job_number, sig))
                #     f.write('VARS\tjob{0}\tfix_source_ratio="{1}"\n'.format(job_number, 'False'))
                #     f.write('VARS\tjob{0}\tsr0="{1}"\n'.format(job_number, 0))
                #     f.write('VARS\tjob{0}\tsr1="{1}"\n'.format(job_number, 0))
                #     f.write('VARS\tjob{0}\tsr2="{1}"\n'.format(job_number, 0))
                #     f.write('VARS\tjob{0}\tfix_scale="{1}"\n'.format(job_number, 'False'))
                #     f.write('VARS\tjob{0}\tscale="{1}"\n'.format(job_number, 0))
                #     f.write('VARS\tjob{0}\tdimension="{1}"\n'.format(job_number, dim))
                #     f.write('VARS\tjob{0}\tenergy="{1}"\n'.format(job_number, en))
                #     f.write('VARS\tjob{0}\tflat_llh="{1}"\n'.format(job_number, flat))
                #     f.write('VARS\tjob{0}\tburnin="{1}"\n'.format(job_number, burnin))
                #     f.write('VARS\tjob{0}\tnwalkers="{1}"\n'.format(job_number, nwalkers))
                #     f.write('VARS\tjob{0}\tnsteps="{1}"\n'.format(job_number, nsteps))
                #     f.write('VARS\tjob{0}\toutfile="{1}"\n'.format(job_number, outchains))
                #     f.write('VARS\tjob{0}\tfix_mixing="{1}"\n'.format(job_number, 'False'))
                #     f.write('VARS\tjob{0}\tno_bsm="{1}"\n'.format(job_number, no_bsm))
                #     job_number += 1
