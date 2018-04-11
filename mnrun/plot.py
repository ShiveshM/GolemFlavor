import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import rc

rc('text', usetex=False)
rc('font', **{'family':'serif', 'serif':['Computer Modern'], 'size':18})

arr = np.load('fr_evidence_test_050_050_000_0100_sfr_033_066_000_DIM3_single_scale.npy')

print arr
ev = arr.T[1]
null = ev[0]
print null
re_z = -(ev - null)
print re_z

plt.plot(arr.T[0], re_z)
plt.xlabel(r'${\rm log}_{10} \Lambda$')
plt.ylabel(r'Bayes Factor')
# plt.axhline(np.power(10, 1.5), color='grey', linestyle='-', alpha=0.4)
plt.savefig('./test.png', bbox_inches='tight', dpi=150)
