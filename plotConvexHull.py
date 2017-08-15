import numpy as np
from scipy.spatial import ConvexHull
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import h5py
from glob import glob

from CSM import init
from CSM import loadData

data = init()
itime = loadData(data, glob('*h5')[0])

points = data['sheep'][itime]

hull = ConvexHull(points)

file = open('./volume','w')
file.write(str(hull.volume))
file.close()

plt.clf()
plt.plot(points[:,0], points[:,1], 'o')
plt.plot(points[hull.vertices,0], points[hull.vertices,1], 'r--', lw=2)
plt.plot(points[hull.vertices[0],0], points[hull.vertices[0],1], 'ro')
plt.axes().set_aspect('equal')
plt.savefig('plots/convexHull.png')
