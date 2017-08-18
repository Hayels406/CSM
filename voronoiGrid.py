import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.spatial import Voronoi, voronoi_plot_2d

coords = np.array([[x,y] for x in np.linspace(-1., 1., 20) for y in np.linspace(-1., 1., 20)])
#plt.scatter(coords[:,0], coords[:,1])
#plt.show()

vor = Voronoi(coords)
voronoi_plot_2d(vor)
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
ax = plt.subplot(111)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(16)

plt.axes().set_aspect('equal')

plt.axis([-1.5, 1.5, -1.5, 1.5])


vertexIndex = np.array(vor.ridge_vertices)[np.where((np.array(vor.ridge_vertices)>=0).sum(axis = 1) == 2)[0]]

dist_v = []

for pair in vertexIndex:
    dist_v.append(np.linalg.norm(vor.vertices[pair[0]]-vor.vertices[pair[1]]))

dist_v = np.array(dist_v)[np.array(dist_v) < .5]

print np.mean(dist_v)

print len(np.where(np.array(vor.ridge_vertices) < 0)[0])

plt.title('Length = ' + str(np.mean(dist_v)) + ', inf = ' + str(len(np.where(np.array(vor.ridge_vertices) < 0)[0])))
plt.show()
