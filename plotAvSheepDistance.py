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
minD = np.array([])
maxD = np.array([])
meanD = np.array([])
t = np.array([])
for dataName in files:
    data = init()
    itime = loadData(data,dataName)

    if dataName == files[0]:
        t = np.append(t, data['t'][0:itime])
        minD = np.append(minD, data['dist_ij']['min'][0:itime])
        maxD = np.append(maxD, data['dist_ij']['max'][0:itime])
        meanD = np.append(meanD, data['dist_ij']['mean'][0:itime])

    else:
        t = np.append(t, data['t'][4:itime])
        minD = np.append(minD, data['dist_ij']['min'][4:itime])
        maxD = np.append(maxD, data['dist_ij']['max'][4:itime])
        meanD = np.append(meanD, data['dist_ij']['mean'][4:itime])

print 'min'
minPlot = plt.figure()
plt.plot(t, minD)
minPlot.show()

print 'max'
maxPlot = plt.figure()
plt.plot(t, maxD)
maxPlot.show()

print 'mean'
meanPlot = plt.figure()
plt.plot(t, meanD)
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

#new.show()
#print data
