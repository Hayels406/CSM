import numpy as np
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys
import h5py
import os
from glob import glob
sys.dont_write_bytecode = True
sys.path.insert(0,os.getcwd())
from params import *


from CSM import init
from CSM import loadData

files = sorted(glob('*.h5'))
t = np.array([])
av = np.array([])
for dataName in files:
    data = init()
    itime = loadData(data,dataName)

    if dataName == files[0]:
        data['t'] = data['t'][0:itime]
        data['sheepAcc'] = data['sheepAcc'][0:itime]
        acc = np.sqrt((data['sheepAcc']**2).sum(axis = 2))[:,0:20]
        maxAcc = np.max(np.sqrt((data['sheepAcc']**2).sum(axis = 2)), axis = 1).tolist()
        tf = len(maxAcc)
        av = np.append(av, np.sqrt((data['sheepAcc']**2).sum(axis = 2)).mean(axis = 1))
    else:
        data['t'] = data['t'][4:itime]
        data['sheepAcc'] = data['sheepAcc'][4:itime]
        acc = np.append(acc, np.sqrt((data['sheepAcc']**2).sum(axis = 2))[:,0:20], axis = 0)
        maxAcc = maxAcc + np.max(np.sqrt((data['sheepAcc']**2).sum(axis = 2)), axis = 1).tolist()
        av = np.append(av, np.sqrt((data['sheepAcc']**2).sum(axis = 2)).mean(axis = 1))
    t = np.append(t, data['t'], axis = 0)

plt.plot(t[10:], acc[10:, :], alpha = 0.2)
plt.plot(t[10:], av[10:,])
plt.plot(t[10:], maxAcc[10:], color = 'r', ls = '--')
plt.yscale('log')
plt.ylabel('log(Acceleration)')
plt.xlabel('Time')
plt.savefig('plots/sheepAcceleration.png')
