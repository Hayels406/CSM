import numpy as np
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys
import h5py
import os
from glob import glob
import matplotlib.mlab as mlab
sys.dont_write_bytecode = True
sys.path.insert(0,os.getcwd())
from params import *


from CSM import init
from CSM import loadData

files = sorted(glob('*.h5'))
acc = np.array([])
for dataName in files:
    data = init()
    itime = loadData(data,dataName)

    if dataName == files[0]:
        data['sheepAcc'] = data['sheepAcc'][200: itime]
        acc = np.append(acc, np.sqrt((data['sheepAcc']**2).sum(axis = 2)))
    else:
        data['sheepAcc'] = data['sheepAcc'][4: itime]
        acc = np.append(acc, np.sqrt((data['sheepAcc']**2).sum(axis = 2)))

plt.figure()
num_bins = 50
mu = acc.mean()
sigma = np.std(acc)
plt.hist(acc, num_bins, normed=1, facecolor='green', alpha=0.5)
n, bins, patches = plt.hist(acc, num_bins, normed=1, facecolor='green', alpha=0.5)
line_best_fit = mlab.normpdf(bins, mu, sigma)
plt.plot(bins, line_best_fit, 'r--')
plt.yscale('log')
plt.savefig('plots/histSheepAcceleration.png')
