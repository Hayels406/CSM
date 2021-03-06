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

area = []

def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

j = 0
for group in sorted(glob('group*[0-9]*')):
	files = glob(group + '/data*.h5')
	print group

	alive = []
	time = []
	plotPeriod = 0.1

	for dFile in files:
		value = dFile[dFile.rfind('-')+1:dFile.rfind('.')]

		data = dict()
		h5f = h5py.File(dFile,'r')
		itime = np.copy(h5f['itime'])[0]
		data['alive'] = np.copy(h5f['alive'])
		data['t'] = np.copy(h5f['t'])

		lastPlot = 0
		for i in range(itime):
			if data['t'][i] - lastPlot > plotPeriod:
				time.append(data['t'][i])
				alive.append(data['alive'][i].sum())
				lastPlot = data['t'][i]


	alive = np.array(alive).reshape(len(alive), 1)
	time = np.array(time).reshape(len(time), 1)

	data = np.append(time, alive, axis = 1)
	data = data.tolist()
	data2 = sorted(data, key=lambda x : x[0])
	data2 = np.array(data2)

	if np.shape(data2)[0] > 0:
		y_av = movingaverage(data2[:,1], 75)
		plt.plot(data2[:,0][100:-50], y_av[100:-50], label = group, color = colors[j])

	trap = np.trapz(y_av[100:-50], x = data2[:,0][100:-50])
	area += [[int(group[5:]), trap]]
	j +=1
plt.colorbar(CS3)
plt.xlabel('Time', fontsize = 18)
plt.ylabel('$N(t)$', fontsize = 18)
plt.savefig('./groupPredation')


np.save('./area.npy', area)
area = np.array(area)
plt.plot(area[:,0], area[:,1], lw = 2)
y_av = movingaverage(area[:,1], 3)
plt.plot(area[:,0][5:-5], y_av[5:-5], lw = 2)
plt.xlabel('Group Size', fontsize = 18)
plt.ylabel('Area', fontsize = 18)
plt.savefig('./groupPredationArea.png')
