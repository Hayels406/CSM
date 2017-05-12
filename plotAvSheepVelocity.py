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
data['sheepVel'] = data['sheepVel'][0:itime]
data['t'] = data['t'][0:itime]
avSpeed = np.sqrt((data['sheepVel']**2).sum(axis = 2)).mean(axis = 1)
plt.plot(data['t'], avSpeed)
plt.show()



#print data
