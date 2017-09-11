import numpy as np
import math
import os
if os.getcwd().rfind('share'):
	import matplotlib as mpl
	mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rc
import sys
import h5py
from glob import glob
import matplotlib.mlab as mlab
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.stats import gamma
sys.dont_write_bytecode = True
sys.path.insert(0,os.getcwd())
from defaultParams import *
from params import *
from fixedParams import *

from CSM import init
from CSM import loadData

dataName = glob('*.h5')[0]

data = init()
itime = loadData(data,dataName)

vor = Voronoi(data['sheep'][itime])
vertexIndex = np.array(vor.ridge_vertices)[np.where((np.array(vor.ridge_vertices)>=0).sum(axis = 1) == 2)[0]]

dist_v = []

for pair in vertexIndex:
    dist_v.append(np.linalg.norm(vor.vertices[pair[0]]-vor.vertices[pair[1]]))

dist_v = np.array(dist_v)[np.array(dist_v) < .5]

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
n, bins, patches = plt.hist(dist_v, bins = 100)
plt.xlabel('Distance between Voronoi vertices', fontsize = 18)
plt.ylabel('Frequency', fontsize = 18)
ax = plt.subplot(111)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(16)
plt.ylim(ymax = 100)
plt.savefig('plots/voronoiDistance.png')
