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
pwd = os.getcwd()
k = int(pwd[pwd.rfind('pred')+4:pwd.rfind('/')])

burnIn = 495

def func(x, alpha):
	t_m = (np.max(time[alive > burnIn]) + np.min(time[alive < burnIn]))/2.
	return NP*np.exp(-alpha*(x-t_m))

def func2(x, beta):
	if np.min(time[alive < k]) == 0.0:
		t_m = (np.max(time[alive > burnIn]) + np.min(time[alive < burnIn]))/2.
	else:
		t_m = (np.max(time[alive > k]) + np.min(time[alive < k]))/2.

	return k*(np.exp(-beta*(x-t_m)))

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

if len(time[alive < k]) == 0:
	t = data['t'][itime]
else:
	t = np.min(time[alive < k])
if t > 0.:
	popt, pcov = curve_fit(func, time[alive > k], alive[alive > k])
	a = popt[0]
	fit = func(np.linspace(0, 100, 1000), a)
	plt.plot(np.linspace(0, 100, 1000)[fit > k], fit[fit > k], label = 'Fit')
else:
	a = 0

if len(time[alive < k]) != 0:
	popt, pcov = curve_fit(func2, time[alive < k], alive[alive < k], p0=[0.01])

	b = popt[0]
	fit = func2(np.linspace(t, 100, 1000), b)
	plt.plot(np.linspace(t, 100, 1000)[fit < k], fit[fit < k], label = 'Fit')
else:
	b = 0

plt.ylim(0,NP +5)
plt.xlim(0,100)
plt.axhline(NP, color = 'red', ls = '--')
plt.axhline(k, color = 'green', ls = ':')


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
plt.close()

file = open('./output/alpha','w')
file.write(str(a))
file.close()
