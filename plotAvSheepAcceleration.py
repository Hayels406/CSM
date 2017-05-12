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
data['sheepAcc'] = data['sheepAcc'][0:itime]
data['t'] = data['t'][0:itime]
av = np.sqrt((data['sheepAcc']**2).sum(axis = 2)).mean(axis = 1)
plt.plot(data['t'], av)
plt.show()



#print data
