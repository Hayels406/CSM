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

for dfile in dataFile:
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
		avV = avV + [np.mean(data['sheepVel'][itime][data['alive'][itime]])]]
		avA = avA + [np.mean(data['sheepAcc'][itime][data['alive'][itime]])]]

	average = average + [np.array(avV).mean(), np.array(avA).mean(), data['dogVel'][:itime].mean(), data['dogAcc'][:itime].mean()]
