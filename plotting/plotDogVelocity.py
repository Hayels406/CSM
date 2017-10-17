import numpy as np
import math
import os
if os.getcwd().rfind('share') > 0:
	topsy = True
	import matplotlib as mpl
	mpl.use('Agg')
else:
	topsy = False
	from matplotlib import rc
import matplotlib.pyplot as plt
import sys
import h5py
from glob import glob
sys.dont_write_bytecode = True
sys.path.insert(0,os.getcwd())
from params import *


from CSM import init
from CSM import loadData

files = sorted(glob('*-1.h5'))
t = np.array([])
av = np.array([])
for dataName in files:
    data = init()
    itime = loadData(data,dataName)

    if dataName == files[0]:
        data['t'] = data['t'][0:itime]
        data['dogVel'] = data['dogVel'][0:itime]
        Vel = np.sqrt((data['dogVel']**2).sum(axis = 1))
        tf = len(Vel)
    else:
        data['t'] = data['t'][4:itime]
        data['dogVel'] = data['dogVel'][4:itime]
        Vel = np.append(Vel, np.sqrt((data['dogVel']**2).sum(axis = 1)), axis = 0)
    t = np.append(t, data['t'], axis = 0)
if topsy == False:
	rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
	rc('text', usetex=True)
plt.plot(t[10:], Vel[10:])
plt.ylim(ymax = np.max(Vel[200:]))
plt.ylabel('Velocity')
plt.xlabel('Time')
plt.savefig('plots/dogVelocity.png')
