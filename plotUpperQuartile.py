import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import rc
import sys
import h5py
import os
from glob import glob
from scipy.spatial import Voronoi, voronoi_plot_2d
sys.dont_write_bytecode = True
sys.path.insert(0,os.getcwd())

wallTypes = ['circular', 'square']

for wall in wallTypes:
    dataFiles = glob(wall+'/*/*h5')

    groupSize = []
    upperQ = []

    if wall == 'circular':
        lab = 'Circular'
    else:
        lab = 'Square'

    for dFile in dataFiles:
        print dFile
        groupSize.append(int(dFile[dFile.rfind('group')+5:dFile.rfind('/')]))

        data = dict()
        h5f = h5py.File(dFile,'r')
    	itime = np.copy(h5f['itime'])[0]
        data['sheep'] = np.copy(h5f['sheep'])

        vor = Voronoi(data['sheep'][itime])
        vertexIndex = np.array(vor.ridge_vertices)[np.where((np.array(vor.ridge_vertices)>=0).sum(axis = 1) == 2)[0]]

        dist_v = []

        for pair in vertexIndex:
            dist_v.append(np.linalg.norm(vor.vertices[pair[0]]-vor.vertices[pair[1]]))

        dist_v = np.array(dist_v)[np.array(dist_v) < .5]

        upperQ.append(np.percentile(dist_v, 75))

    plt.plot(groupSize, upperQ, lw = 2, label = lab)
    plt.scatter(groupSize, upperQ, color = 'Black')
    if wall == 'circular':
        groupSizeC = groupSize
        upperQC = upperQ
    else:
        groupSizeS = groupSize
        upperQS = upperQ

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
plt.xlabel('Group size', fontsize = 18)
plt.ylabel('Upper Quartile of Voronoi Tesselation Lengths', fontsize = 18)
plt.xlim(-5,405)
plt.legend()
ax = plt.subplot(111)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(16)
plt.savefig('./upperQ.png')
np.savetxt('./upperQ.out', (groupSizeC,upperQC,groupSizeS,upperQS))
