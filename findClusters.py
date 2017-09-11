import numpy as np
import h5py
from glob import glob
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
import os
if os.getcwd().rfind('share'):
	import matplotlib as mpl
	mpl.use('Agg')
import matplotlib.pyplot as plt

dFile = glob('*h5')[0]
av_score = []
pwd = os.getcwd()

data = dict()
h5f = h5py.File(dFile,'r')
itime = np.copy(h5f['itime'])[0]
data['sheep'] = np.copy(h5f['sheep'])
NP = 500.
groupSize = int(pwd[pwd.rfind('/group')+len('/group'):])

if int(NP/groupSize) != 1:
    for n_clusters in np.array(range(int(NP/groupSize) - 1))+2:
        clusterer = KMeans(n_clusters=n_clusters, random_state=10)
        cluster_labels = clusterer.fit_predict(data['sheep'][itime])

        print "For n_clusters =", n_clusters, "the average silhouette_score is :", silhouette_score(data['sheep'][itime], cluster_labels)

        av_score = av_score + [silhouette_score(data['sheep'][itime], cluster_labels)]

    if np.max(av_score) < 0.5:
        num_clusters = 1
        plt.scatter(data['sheep'][itime][:,0], data['sheep'][itime][:,1])
        plt.axes().set_aspect('equal')
        plt.savefig('plots/clusters.png')
        file = open('./output/clusters','w')
        file.write(str(1))
        file.close()
    else:
        num_clusters = np.where(av_score == np.max(av_score))[0][0]+2

        print 'There are', num_clusters, 'clusters'
        clusterer = KMeans(n_clusters=num_clusters, random_state=10)
        cluster_labels = clusterer.fit_predict(data['sheep'][itime])

        for i in range(num_clusters):
            plt.scatter(data['sheep'][itime][cluster_labels == i][:,0], data['sheep'][itime][cluster_labels == i][:,1])
            plt.axes().set_aspect('equal')
        plt.savefig('plots/clusters.png')
        file = open('./output/clusters','w')
        file.write(str(num_clusters))
        file.close()

else:
    plt.scatter(data['sheep'][itime][:,0], data['sheep'][itime][:,1])
    plt.axes().set_aspect('equal')
    plt.savefig('plots/clusters.png')
    file = open('./output/clusters','w')
    file.write(str(1))
    file.close()
