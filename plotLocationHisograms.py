import numpy as np
import math
import matplotlib.pyplot as plt
import sys
import h5py
import os
from glob import glob
sys.dont_write_bytecode = True
sys.path.insert(0,os.getcwd())
from defaultParams import *
from fixedParams import *
from params import *

from CSM import init
from CSM import loadData

overTime = 'Off'
meanOverTime = 'Off'
finalTime = 'On'

files = sorted(glob('*.h5'))
plt.figure()
for dataName in files:
    data = init()
    itime = loadData(data,dataName)
    lastPlot = data['t'][0]
    tstep = 0
    if overTime == 'On':
        for t in data['t'][:itime]:
            if t-lastPlot > plotPeriod:
                print t
                plt.clf()
                plt.hist(data['sheep'][:itime,:,0].mean(axis = 0), bins = 150, range = (-1,1), label = "Mean over time")
                #plt.hist(data['sheep'][:itime,:,0], bins = 100, range = (-1,1), label = "Cumulative")
                #plt.hist(data['sheep'][tstep,:,0], bins = 150, range = (-1,1),alpha = 0.9, label = "At t = " + str(data['t'][tstep]))
                lastPlot = data['t'][tstep]
                plt.legend()
                plt.ylim(ymax = 9)
                plt.xlabel('Mean position in $x$')
                plt.ylabel('Frequency')
                plt.pause(0.05)
            tstep += 1

    elif meanOverTime == 'On':
        plt.hist(data['sheep'][:itime,:,0].mean(axis = 0), bins = 150, range = (-1,1), label = "Mean over time")
        plt.ylim(ymax = 9)
        plt.xlabel('Mean position in $x$')
        plt.ylabel('Frequency')
        plt.savefig('plots/meanXLocationHistogram.png')

    elif finalTime == 'On':
        plt.hist(data['sheep'][itime,:,0], bins = 150, range = (-1,1))
        plt.ylim(ymax = 9)
        plt.xlabel('Final position in $x$')
        plt.ylabel('Frequency')
        plt.savefig('plots/finalXLocationHistogram.png')
