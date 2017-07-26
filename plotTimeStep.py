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
from defaultParams import *
from params import *


from CSM import init
from CSM import loadData

files = sorted(glob('*.h5'))
t = np.array([])
q = np.array([])
for dataName in files:
    data = init()
    itime = loadData(data,dataName)

    if dataName == files[0]:
        data['t'] = data['t'][3:itime]
        data['q'] = data['q'][0:itime]
    else:
        data['t'] = data['t'][4:itime]
        data['q'] = data['q'][1:itime]
    t = np.append(t, data['t'], axis = 0)
    q = np.append(q, data['q'], axis = 0)

plt.plot(t[10:], dt*q[10:])
plt.ylabel('Time Step Size')
plt.xlabel('Time')
plt.savefig('plots/timeStep.png')
