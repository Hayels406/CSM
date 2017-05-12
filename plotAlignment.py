import numpy as np
import math
import matplotlib.pyplot as plt
import sys
import h5py
import os
from scipy.spatial.distance import cdist
sys.dont_write_bytecode = True
sys.path.insert(0,os.getcwd())
from params import *


from CSM import init
from CSM import loadData

data = init()
itime = loadData(data,fileName)

data['t'] = data['t'][0:itime]
data['sheepVel'] = data['sheepVel'][0:itime]

tstep = 0
align = np.zeros(itime)
for t in data['t']:
    thetas = np.arctan2(data['sheepVel'][tstep][:,0], data['sheepVel'][tstep][:,1])
    math.sqrt(pow(np.cos(thetas).mean(), 2) + pow(np.sin(thetas).mean(), 2))
    align[tstep] = math.sqrt(pow(np.cos(thetas).mean(), 2) + pow(np.sin(thetas).mean(), 2))
    tstep += 1

plt.plot(data['t'], align)
plt.ylim(0, 1)
plt.show()
