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

from params import *

dataFile = glob('data*-*.h5')
average = []

for dFile in dataFile:
	data = dict()
	h5f = h5py.File(dFile,'r')
	itime = np.copy(h5f['itime'])[0]
	data['alive'] = np.copy(h5f['alive'])
	data['sheepAcc'] = np.copy(h5f['sheepAcc'])
	data['sheepVel'] = np.copy(h5f['sheepVel'])
	data['dogAcc'] = np.copy(h5f['dogAcc'])
	data['dogVel'] = np.copy(h5f['dogVel'])

	#sheep
	avV = []
	avA = []
	for i in range(itime):
		avV = avV + [[  np.mean(np.sqrt((data['sheepVel'][i][data['alive'][i]]**2).sum(axis = 1)))   ]]
		avA = avA + [[np.mean(np.sqrt((data['sheepAcc'][i][data['alive'][i]]**2).sum(axis = 1)))]]

	average = average + [[np.array(avV).mean(), np.array(avA).mean(), np.mean(np.sqrt((data['dogVel']**2).sum(axis = 1))),             np.mean(np.sqrt((data['dogAcc']**2).sum(axis = 1)))]]

average = np.array(average)

plt.plot(average[:,0], label = 'Prey Vel')
plt.scatter(range(10), average[:,0])
plt.plot(average[:,1], label = 'Prey Acc')
plt.scatter(range(10), average[:,1])
plt.plot(average[:,2], label = 'Pred Vel')
plt.scatter(range(10), average[:,2])
plt.plot(average[:,3], label = 'Pred Acc')
plt.scatter(range(10), average[:,3])
plt.legend(loc = 'best')
plt.savefig('del.png')
