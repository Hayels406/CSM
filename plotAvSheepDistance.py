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
#from CSM import loadData

#def avSpeed(data):
#    (np.sqrt((data['sheepVel']**2).sum(axis = 2))).mean(axis = 1)

data = init()
itime = loadData(data,fileName)
data['t'] = data['t'][0:itime]
data['dist_ij']['min'] = data['dist_ij']['min'][0:itime]
data['dist_ij']['max'] = data['dist_ij']['max'][0:itime]
data['dist_ij']['mean'] = data['dist_ij']['mean'][0:itime]

minPlot = plt.figure()
plt.plot(data['t'], data['dist_ij']['min'])
minPlot.show()

maxPlot = plt.figure()
plt.plot(data['t'], data['dist_ij']['max'])
maxPlot.show()

meanPlot = plt.figure()
plt.plot(data['t'], data['dist_ij']['mean'])
meanPlot.show()

#new  = plt.figure()
#from mpl_toolkits.axes_grid1 import AxesGrid
#gr = AxesGrid(plt.gcf(),111,[3,1], axes_pad=0,)
#gr[0].plot(data['t'], data['dist_ij']['min'])
#gr[1].plot(data['t'], data['dist_ij']['max'])
#gr[1].locator_params(axis = 'y', prune = 'upper')
#gr[2].plot(data['t'], data['dist_ij']['mean'])
#gr[2].locator_params(axis = 'y', prune = 'upper')
#for g in gr:
#    g.set_ylim(0.9, 10)

new.show()
#print data
