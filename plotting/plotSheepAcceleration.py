import numpy as np
import math
import os
if os.getcwd().rfind('share') > 0:
	topsy = True
	import matplotlib as mpl
	mpl.use('Agg')
else:
	topsy = False
	from matplotlib import rc

import matplotlib.pyplot as plt
import sys
import h5py
from glob import glob
sys.dont_write_bytecode = True
sys.path.insert(0,os.getcwd())
from params import *


from CSM import init
from CSM import loadData

files = sorted(glob('*-1.h5'))
t = np.array([])
av = np.array([])
maxAcc = np.array([])
for dataName in files:
	data = init()
	h5f = h5py.File(dFile,'r')
	itime = np.copy(h5f['itime'])[0]
	data['alive'] = np.copy(h5f['alive'])
	data['t'] = np.copy(h5f['t'])
	data['sheepAcc'] = np.copy(h5f['sheepAcc'])

    if dataName == files[0]:
        data['t'] = data['t'][0:itime]
        data['sheepAcc'] = data['sheepAcc'][0:itime]
        acc = np.sqrt((data['sheepAcc']**2).sum(axis = 2))[:,data['alive'][itime]][:,0:20]
        for i in range(itime):
            maxAcc = np.append(maxAcc, np.max(np.sqrt((data['sheepAcc'][i][data['alive'][i]]**2).sum(axis = 1)), axis = 0))
            av = np.append(av, np.sqrt((data['sheepAcc'][i][data['alive'][i]]**2).sum(axis = 1)).mean(axis = 0))
    else:
        data['t'] = data['t'][4:itime]
        data['sheepAcc'] = data['sheepAcc'][4:itime]
        acc = np.append(acc, np.sqrt((data['sheepAcc']**2).sum(axis = 2))[:,data['alive'][itime]][:,0:20], axis = 0)
        for i in range(itime):
            maxAcc = np.append(maxAcc, np.max(np.sqrt((data['sheepAcc'][i][data['alive'][i]]**2).sum(axis = 1)), axis = 0))
            av = np.append(av, np.sqrt((data['sheepAcc'][i][data['alive'][i]]**2).sum(axis = 1)).mean(axis = 0))

    t = np.append(t, data['t'], axis = 0)

plt.plot(t[t<=2000][10:], acc[t<=2000][10:, :], alpha = 0.2)
plt.plot(t[t<=2000][10:], av[t<=2000][10:,])
plt.plot(t[t<=2000][10:], np.array(maxAcc)[t<=2000][10:], color = 'r', ls = '--')
if topsy == False:
	rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
	rc('text', usetex=True)
plt.yscale('log')
plt.ylabel('log(Acceleration)')
plt.xlabel('Time')
ax = plt.subplot(111)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(16)
plt.savefig('plots/sheepAcceleration.png')
plt.close()
