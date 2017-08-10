import numpy as np
import math
import matplotlib.pyplot as plt
import sys
import h5py
import os
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

n, bins, patches = plt.hist(dist_v, bins = 100)
plt.xlabel('Distance between Voronoi vertices')
plt.ylabel('Frequency')
#mu = np.mean(dist_v)
#sd = np.std(dist_v)
#beta = mu/sd**2
#alpha = mu**2/sd**2
#fit_alpha, fit_loc, fit_beta=gamma.fit(dist_v)
#line_best_fit = gamma.pdf(bins, scale = 1./beta, a = alpha)
#plt.plot(bins, line_best_fit)
plt.ylim(ymax = 100)
plt.title('mean = ' + str(np.mean(dist_v)) + '   Upper Quartile = ' + str(np.percentile(dist_v, 75)))
plt.savefig('plots/voronoiDistance.png')
