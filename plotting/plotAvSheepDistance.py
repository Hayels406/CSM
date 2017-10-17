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
plt.plot(t, minD, lw = 2)
plt.axhline(mean(minD), color = 'r', lw = 2)
minPlot.savefig('minDistance.png')

print 'max'
maxPlot = plt.figure()
plt.plot(t, maxD, lw = 2)
plt.axhline(mean(maxD), color = 'r', lw = 2)
maxPlot.savefig('minDistance.png')

print 'mean'
meanPlot = plt.figure()
plt.plot(t, meanD, lw = 2)
plt.axhline(mean(meanD), color = 'r', lw = 2)
meanPlot.savefig('meanDisance.png')
