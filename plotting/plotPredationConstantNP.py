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
from matplotlib import cm
from scipy.optimize import curve_fit
import sys
import h5py
from glob import glob

pwd = os.getcwd()
k = int(pwd[pwd.rfind('pred')+4:])

number_of_lines = len(glob('group*[0-9]*'))
cm_subsection = np.linspace(0., 1., number_of_lines)
colors = [ cm.magma(x) for x in cm_subsection]

Z = [[0,0],[0,0]]
levels = range(5,500+5,5)
CS3 = plt.contourf(Z, levels, cmap='magma')
plt.clf()

a = []

def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

def linear_fit(x, a, b):
	return a*x + b

j = 0
for group in sorted(glob('group*[0-9]*')):
	files = glob(group + '/data*.h5')
	print group

	deathCount = []
	time = []
	plotPeriod = 0.1

	for dFile in files:
		value = dFile[dFile.rfind('-')+1:dFile.rfind('.')]

		data = dict()
		h5f = h5py.File(dFile,'r')
		itime = np.copy(h5f['itime'])[0]
		data['deathCount'] = np.copy(h5f['deathCount'])
		data['t'] = np.copy(h5f['t'])

		lastPlot = 0
		for i in range(itime):
			if data['t'][i] - lastPlot > plotPeriod:
				time.append(data['t'][i])
				deathCount.append(data['deathCount'][i].sum())
				lastPlot = data['t'][i]


	deathCount = np.array(deathCount).reshape(len(deathCount), 1)
	time = np.array(time).reshape(len(time), 1)

	data = np.append(time, deathCount, axis = 1)
	data = data.tolist()
	data2 = sorted(data, key=lambda x : x[0])
	data2 = np.array(data2)

	if np.shape(data2)[0] > 0:
		y_av = movingaverage(data2[:,1], 75)
		popt, pcov = curve_fit(linear_fit, data2[:,0][100:-50], y_av[100:-50])
		grad = popt[0]
		plt.plot(data2[:,0][100:-50], y_av[100:-50], label = group, color = colors[j])

	a += [[int(group[5:]), grad]]
	j +=1
plt.colorbar(CS3)
plt.xlabel('Time', fontsize = 18)
plt.ylabel('$N(t)$', fontsize = 18)
plt.savefig('./predationConstantNP.png')

plt.close()
a = np.array(a)
plt.plot(a[:,0], a[:,1], lw = 2, label = 'Gradient')
y_av = movingaverage(a[:,1], 5)
plt.plot(a[:,0][5:-5], y_av[5:-5], lw = 2, label = 'Moving Average')
plt.legend(loc = 'best')
plt.xlabel('Group Size', fontsize = 18)
plt.ylabel('Gradient', fontsize = 18)
plt.savefig('./predationConstantNPGradient.png')
