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


files = sorted(glob('*-1.h5'))
t = np.array([])
for dataName in files:
	data = dict()
	h5f = h5py.File(dataName,'r')
	itime = np.copy(h5f['itime'])[0]
	data['alive'] = np.copy(h5f['alive'])
	data['t'] = np.copy(h5f['t'])
	data['dogAcc'] = np.copy(h5f['dogAcc'])

	if dataName == files[0]:
		data['t'] = data['t'][0:itime]
		data['dogAcc'] = data['dogAcc'][0:itime]
		Acc = np.sqrt((data['dogAcc']**2).sum(axis = 1))
		tf = len(Acc)
	else:
		data['t'] = data['t'][4:itime]
		data['dogAcc'] = data['dogAcc'][4:itime]
		Acc = np.append(Acc, np.sqrt((data['dogAcc']**2).sum(axis = 1)), axis = 0)
	t = np.append(t, data['t'], axis = 0)
if topsy == False:
	rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
	rc('text', usetex=True)
plt.plot(t[10:], Acc[10:])
plt.ylim(ymax = np.max(Acc[200:]))
plt.ylabel('Acceleration')
plt.xlabel('Time')
plt.savefig('plots/dogAcceleration.png')
