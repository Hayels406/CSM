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
from scipy.optimize import curve_fit
import sys
import h5py
from glob import glob

if topsy == False:
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

filesSquare = glob('./pred*/group*/output/beta')
filesSquare.sort()

k = np.array([])
Ni = np.array([])
beta = np.array([])
for fileLoc in filesSquare:
	k = np.append(k, int(fileLoc[fileLoc.rfind('pred')+4:fileLoc.rfind('/group')]))
	Ni = np.append(Ni, int(fileLoc[fileLoc.rfind('group')+5:fileLoc.rfind('output')-1]))
	file = open(fileLoc,'r')
	betaValue = float(file.read())
	file.close()
	beta = np.append(beta, betaValue)

k = np.array(map(int, k))
Ni = np.array(map(int, Ni))
beta = beta.reshape(len(Ni)/((Ni == 500).sum()), (Ni == 500).sum())

plt.pcolor(beta)
plt.colorbar()
plt.xlabel('Pred Neighbours', fontsize = 18)
plt.ylabel('Prey Neighbours', fontsize = 18)
ax = plt.gca()
ax.set_xticks(np.array(range(len(list(set(k))))) + 0.5)
ax.set_xticklabels(sorted(list(set(k))))
ax.set_yticks(np.array(range(len(list(set(Ni))))) + 0.5)
ax.set_yticklabels(sorted(list(set(Ni))))
plt.savefig('./pcolorBeta.png')

plt.close()
for i in range(len(set(k))):
	kValue = np.array(sorted(list(set(k))))[i]
	plt.plot(Ni[k == kValue], beta[:, i])
	plt.xlabel('Prey Neighbours', fontsize = 18)
	plt.ylabel(r'$\beta$', fontsize = 18)
	plt.title(r'$k = '+str(kValue)+'$', fontsize = 20)
	plt.savefig('./pred'+str(kValue).zfill(3)+'/beta.png')
	plt.close()
