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
maxVel = np.array([])
for dataName in files:
	data = dict()
	h5f = h5py.File(dataName,'r')
	itime = np.copy(h5f['itime'])[0]
	data['alive'] = np.copy(h5f['alive'])
	data['t'] = np.copy(h5f['t'])
	data['sheepVel'] = np.copy(h5f['sheepVel'])

	if dataName == files[0]:
		data['t'] = data['t'][0:itime]
		data['sheepVel'] = data['sheepVel'][0:itime]
		Vel = np.sqrt((data['sheepVel']**2).sum(axis = 2))[:,data['alive'][itime]][:,0:20]
		for i in range(itime):
			maxVel = np.append(maxVel, np.max(np.sqrt((data['sheepVel'][i][data['alive'][i]]**2).sum(axis = 1)), axis = 0))
			av = np.append(av, np.sqrt((data['sheepVel'][i][data['alive'][i]]**2).sum(axis = 1)).mean(axis = 0))
	else:
		data['t'] = data['t'][4:itime]
		data['sheepVel'] = data['sheepVel'][4:itime]
		Vel = np.append(Vel, np.sqrt((data['sheepVel']**2).sum(axis = 1))[data['alive'][itime]][:,0:20], axis = 0)
		for i in range(itime):
			maxVel = np.append(maxVel, np.max(np.sqrt((data['sheepVel'][i][data['alive'][i]]**2).sum(axis = 1)), axis = 0))
			av = np.append(av, np.sqrt((data['sheepVel'][i][data['alive'][i]]**2).sum(axis = 1)).mean(axis = 0))
	t = np.append(t, data['t'], axis = 0)
if topsy == False:
	rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
	rc('text', usetex=True)
plt.plot(t[10:], Vel[10:, :], alpha = 0.2)
plt.plot(t[10:], av[10:,])
plt.plot(t[10:], maxVel[10:], color = 'r', ls ='--')
plt.ylim(ymax = np.max(maxVel[200:]))
plt.ylabel('Velocity')
plt.xlabel('Time')
plt.savefig('plots/sheepVelocity.png')
