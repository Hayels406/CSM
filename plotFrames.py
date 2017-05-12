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
from CSM import initPlot
from CSM import plotDataPositions

data = init()
itime = loadData(data,fileName)

data['t'] = data['t'][0:itime]
data['sheep'] = data['sheep'][0:itime]
data['sheepVel'] = data['sheepVel'][0:itime]
data['dog'] = data['dog'][0:itime]
data['dogVel'] = data['dogVel'][0:itime]

tstep = 0
lastPlot = data['t'][tstep]
dogQuiver, sheepQuiver = initPlot(data)
for t in data['t']:
    if t-lastPlot > plotPeriod:
		print t
		plotDataPositions(data, tstep, dogQuiver, sheepQuiver)
		lastPlot = data['t'][tstep]
    tstep += 1
