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
av = np.array([])
t = np.array([])
for dataName in files:
    data = init()
    itime = loadData(data,dataName)

    if dataName == files[0]:
        data['t'] = data['t'][0:itime]
        data['sheepAcc'] = data['sheepAcc'][0:itime]
        av = np.append(av, np.sqrt((data['sheepAcc']**2).sum(axis = 2)).mean(axis = 1))
    else:
        data['t'] = data['t'][4:itime]
        data['sheepAcc'] = data['sheepAcc'][4:itime]
        av = np.append(av, np.sqrt((data['sheepAcc']**2).sum(axis = 2)).mean(axis = 1))
    t = np.append(t, data['t'])


plt.plot(t, av)
plt.show()
