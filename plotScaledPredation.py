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

dataFile = glob('group*/*.h5')
save = []
for dFile in dataFile:
    data = dict()
    h5f = h5py.File(dFile,'r')
    itime = np.copy(h5f['itime'])[0]
    data['alive'] = np.copy(h5f['alive'])
    data['t'] = np.copy(h5f['t'])
    t = ['t_' + str(int(dFile[dFile.rfind('p')+1:dFile.rfind('/')]))]
    alive = ['alive_' + str(int(dFile[dFile.rfind('p')+1:dFile.rfind('/')]))]
    plt.plot(data['t'][:itime], (1./data['alive'][0].sum())*data['alive'][:itime].sum(axis = 1), lw = 2, label = 'Group size = ' + str(int(dFile[dFile.rfind('p')+1:dFile.rfind('/')])))
    t = t + map(str, data['t'][:itime])
    alive = alive + map(str, (1./data['alive'][0].sum())*data['alive'][:itime].sum(axis = 1)
)
    save = save + [t] + [alive]

plt.ylim(0,1.01)
plt.axhline(1, color = 'red', ls = '--')
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
plt.xlabel('Time', fontsize = 18)
plt.ylabel('$N(t)/N(0)$', fontsize = 18)
ax = plt.subplot(111)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(16)
plt.legend(loc = 'upper right')
plt.savefig('./scaledPredation.png')
maxLen = 0
for i in range(np.shape(save)[0]):
    if maxLen < len(save[i]):
        maxLen = len(save[i])

for i in range(np.shape(save)[0]):
    while len(save[i]) < maxLen:
        save[i].append('')

np.savetxt('./scaledPredation.csv', np.transpose(np.array(save)), fmt = '%s', delimiter = ',')
