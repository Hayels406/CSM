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
from matplotlib import rc
from scipy.optimize import curve_fit
import sys
import h5py
from glob import glob

dataFile = glob('data*-*.h5')

def func(x, b):
    return NP*np.exp(-b*x)

alive = np.array([])
time = np.array([])

for dFile in dataFile:
	print dFile
	value = dFile[dFile.rfind('-')+1:dFile.rfind('.')]

	data = dict()
	h5f = h5py.File(dFile,'r')
	itime = np.copy(h5f['itime'])[0]
	data['alive'] = np.copy(h5f['alive'])
	data['t'] = np.copy(h5f['t'])

	alive = np.append(alive, data['alive'][:itime].sum(axis = 1))
	time = np.append(time, data['t'][:itime])
	if dFile == dataFile[0]:
		NP = np.shape(data['alive'])[1]
		plt.plot(data['t'][:itime], data['alive'][:itime].sum(axis = 1), lw = 2, color = 'k', alpha = 0.5, label = 'Iterations')
	else:
		plt.plot(data['t'][:itime], data['alive'][:itime].sum(axis = 1), lw = 2, color = 'k', alpha = 0.5)


popt, pcov = curve_fit(func, time, alive)
b = popt[0]

plt.ylim(0,NP +5)
plt.xlim(0,100)
plt.axhline(NP, color = 'red', ls = '--')
plt.plot(np.linspace(0, 100, 1000), func(np.linspace(0, 100, 1000), b), label = 'Fit')
if topsy == False:
	rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
	rc('text', usetex=True)
plt.xlabel('Time', fontsize = 18)
plt.ylabel('$N(t)$')
plt.legend(loc = 'upper right', fontsize = 16)
plt.savefig('./plots/ensemblePredation.png')

file = open('./output/beta','w')
file.write(str(b))
file.close()
