import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import rc
import sys
import h5py
import os
from glob import glob

from CSM import init
from CSM import loadData

data = init()
itime = loadData(data, glob('*h5')[0])

plt.plot(data['t'][:itime], data['alive'][:itime].sum(axis = 1), lw = 2)
plt.ylim(0,405)
plt.axhline(400, color = 'red', ls = '--')
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
plt.xlabel('Time', fontsize = 18)
plt.ylabel('$N(t)$', fontsize = 18)
ax = plt.subplot(111)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(16)
plt.savefig('./plots/predation.png')
np.savetxt('./output/predation.out', [data['t'][:itime], data['alive'][:itime].sum(axis = 1)])
