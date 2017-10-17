import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import rc
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
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

for dataName in files:
    data = init()
    itime = loadData(data,dataName)
    lastPlot = data['t'][0]
    tstep = 0

    plt.hist(data['sheep'][itime,:,0], bins = 150, range = (-1.25,1.25))
    plt.ylim(ymax = 9)
    plt.xlim(-1.25, 1.25)
    plt.xlabel('Final position in $x$', fontsize = 18)
    plt.ylabel('Frequency', fontsize = 18)
    ax = plt.subplot(111)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(16)
    plt.savefig('plots/finalXLocationHistogram.png')

    plt.clf()
    plt.hist(data['sheep'][:itime,:,0].mean(axis = 0), bins = 150, range = (-1.25,1.25), label = "Mean over time")
    plt.ylim(ymax = 9)
    plt.xlim(-1.25, 1.25)
    plt.xlabel('Mean position in $x$', fontsize = 18)
    plt.ylabel('Frequency', fontsize = 18)
    ax = plt.subplot(111)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(16)
    plt.savefig('plots/meanXLocationHistogram.png')

#    plt.clf()
#    for t in data['t'][:itime]:
#        if t-lastPlot > plotPeriod:
#            print t
#            plt.clf()
#            plt.hist(data['sheep'][:itime,:,0].mean(axis = 0), bins = 150, range = (-1,1), label = "Mean over time")
#                #plt.hist(data['sheep'][:itime,:,0], bins = 100, range = (-1,1), label = "Cumulative")
#                #plt.hist(data['sheep'][tstep,:,0], bins = 150, range = (-1,1),alpha = 0.9, label = "At t = " + str(data['t'][tstep]))
#            lastPlot = data['t'][tstep]
#            plt.xlim(-1.5, 1.5)
#            plt.ylim(ymax = 9)
#            plt.xlabel('Mean position in $x$')
#            plt.ylabel('Frequency')
#            plt.pause(0.05)
#        tstep += 1
