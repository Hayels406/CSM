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
for dataName in files:
    data = init()
    itime = loadData(data,dataName)

    if dataName == files[0]:
        data['t'] = data['t'][0:itime]
        data['dogAcc'] = data['dogAcc'][0:itime]
        Acc = np.sqrt((data['dogAcc']**2).sum(axis = 1))
        tf = len(Acc)
    else:
        data['t'] = data['t'][4:itime]
        data['dogAcc'] = data['dogAcc'][4:itime]
        Acc = np.append(Acc, np.sqrt((data['dogAcc']**2).sum(axis = 1)), axis = 0)
    t = np.append(t, data['t'], axis = 0)

plt.plot(t[10:], Acc[10:])
plt.ylim(ymax = np.max(Acc[200:]))
plt.ylabel('Acceleration')
plt.xlabel('Time')
plt.savefig('plots/dogAcceleration.png')
