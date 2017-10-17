import numpy as np
import math
import sys
import h5py
from glob import glob
import os
if os.getcwd().rfind('share'):
	import matplotlib as mpl
	mpl.use('Agg')
import matplotlib.pyplot as plt
sys.dont_write_bytecode = True
sys.path.insert(0,os.getcwd())

from defaultParams import *
from fixedParams import *
from params import *


from CSM import init
from CSM import loadData
from CSM import initPlot
from CSM import plotDataPositions

if os.getcwd().rfind('share') > 0:
	import matplotlib as mpl
	mpl.use('Agg')
	topsy = True
else:
	topsy = False


files = sorted(glob('*-1.h5'))

for dataName in files:
    data = init()
    itime = loadData(data,dataName)

    if dataName == files[0]:
        data['t'] = data['t'][0:itime]
        data['sheep'] = data['sheep'][0:itime]
        data['sheepVel'] = data['sheepVel'][0:itime]
        data['dog'] = data['dog'][0:itime]
        data['dogVel'] = data['dogVel'][0:itime]
        if segments == 'On':
            data['interactingSheep'] = data['interactingSheep'][0:itime]
    else:
        data['t'] = data['t'][4:itime]
        data['sheep'] = data['sheep'][4:itime]
        data['sheepVel'] = data['sheepVel'][4:itime]
        data['dog'] = data['dog'][4:itime]
        data['dogVel'] = data['dogVel'][4:itime]
        if segments == 'On':
            data['interactingSheep'] = data['interactingSheep'][4:itime]

    tstep = 0
    lastPlot = data['t'][tstep]
    if dataName == files[0]:
        if predOff == False:
            dogQuiver, sheepQuiver = initPlot(data, savePlotPng)
        else:
            sheepQuiver = initPlot(data, savePlotPng)
    for t in data['t']:
        if t-lastPlot > plotPeriod:
            print 'read', t
            if predOff == False:
                plotDataPositions(data, tstep, dogQuiver, sheepQuiver, savePlotPng)
            else:
                plotDataPositions(data, tstep, 'off', sheepQuiver, savePlotPng)
    		lastPlot = data['t'][tstep]
        tstep += 1
