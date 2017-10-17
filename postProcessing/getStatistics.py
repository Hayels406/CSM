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
AVERAGE = np.mean(average, axis =0)
file = open('./output/predVel','w')
file.write(str(AVERAGE[2]))
file.close()
file = open('./output/predAcc','w')
file.write(str(AVERAGE[3]))
file.close()
file = open('./output/preyVel','w')
file.write(str(AVERAGE[0]))
file.close()
file = open('./output/preyAcc','w')
file.write(str(AVERAGE[1]))
file.close()
