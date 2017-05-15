import numpy as np
import math
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
align = np.array([])
t = np.array([])
for dataName in files:
    data = init()
    itime = loadData(data,dataName)

    if dataName == files[0]:
        data['t'] = data['t'][0:itime]
        data['sheepVel'] = data['sheepVel'][0:itime]
        loop = 0
        align = np.append(align, np.zeros(itime))
    else:
        data['t'] = data['t'][4:itime]
        data['sheepVel'] = data['sheepVel'][4:itime]
        loop = loop + tstep
        align = np.append(align, np.zeros(np.shape(data['t'])[0]))
    t = np.append(t, data['t'])


    for tstep in range(len(data['t'])):
        thetas = np.arctan2(data['sheepVel'][tstep][:,0], data['sheepVel'][tstep][:,1])
        math.sqrt(pow(np.cos(thetas).mean(), 2) + pow(np.sin(thetas).mean(), 2))
        align[tstep + loop] = math.sqrt(pow(np.cos(thetas).mean(), 2) + pow(np.sin(thetas).mean(), 2))

plt.plot(t, align, lw = 2)
plt.ylim(-0.005, 1.05)
plt.axhline(1.0, color = 'r', lw = 2)
plt.show()

if sys.argv[1] == 'Save':
    plt.savefig('alignment.png')
