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

filesSquare = glob('Square/pred200/group*/output/alpha')
filesSquare.sort()

filesCircular = glob('Circular/pred200/group*/output/alpha')
filesCircular.sort()

k = np.array([])
Ni = np.array([])
alpha = np.array([])
for fileLoc in filesSquare:
	k = np.append(k, int(fileLoc[fileLoc.rfind('pred')+4:fileLoc.rfind('/group')]))
	Ni = np.append(Ni, int(fileLoc[fileLoc.rfind('group')+5:fileLoc.rfind('output')-1]))
	file = open(fileLoc,'r')
	betaValue = float(file.read())
	file.close()
	alpha = np.append(alpha, betaValue)

k = np.array(map(int, k))
Ni = np.array(map(int, Ni))
alpha = alpha.reshape(len(Ni)/((Ni == 500).sum()), (Ni == 500).sum())

plt.pcolor(alpha)
plt.colorbar()
plt.xlabel('Pred Neighbours', fontsize = 18)
plt.ylabel('Prey Neighbours', fontsize = 18)
ax = plt.gca()
ax.set_xticks(np.array(range(len(list(set(k))))) + 0.5)
ax.set_xticklabels(sorted(list(set(k))))
ax.set_yticks(np.array(range(len(list(set(Ni))))) + 0.5)
ax.set_yticklabels(sorted(list(set(Ni))))
plt.savefig('./Square/pcolorAlpha.png')

plt.close()
for i in range(len(set(k))):
	kValue = np.array(sorted(list(set(k))))[i]
	plt.plot(Ni[k == kValue], alpha[:, i])
	plt.xlabel('Prey Neighbours', fontsize = 18)
	plt.ylabel(r'$\alpha$', fontsize = 18)
	plt.title(r'$k = '+str(kValue)+'$', fontsize = 20)
	plt.savefig('./Square/pred'+str(kValue).zfill(3)+'/alpha.png')
	plt.close()


k = np.array([])
Ni = np.array([])
alpha = np.array([])
for fileLoc in filesCircular:
	fileLoc[fileLoc.rfind('/group')]
	k = np.append(k, int(fileLoc[fileLoc.rfind('pred')+4:fileLoc.rfind('/group')]))
	Ni = np.append(Ni, int(fileLoc[fileLoc.rfind('group')+5:fileLoc.rfind('output')-1]))
	file = open(fileLoc,'r')
	betaValue = float(file.read())
	file.close()
	alpha = np.append(alpha, betaValue)

k = np.array(map(int, k))
Ni = np.array(map(int, Ni))
alpha = alpha.reshape(len(Ni)/((Ni == 500.).sum()), (Ni == 500.).sum())


plt.pcolor(alpha)
plt.colorbar()
plt.xlabel('Pred Neighbours')
plt.ylabel('Prey Neighbours')
ax = plt.gca()
ax.set_xticks(np.array(range(len(list(set(k))))) + 0.5)
ax.set_xticklabels(sorted(list(set(k))))
ax.set_yticks(np.array(range(len(list(set(Ni))))) + 0.5)
ax.set_yticklabels(sorted(list(set(Ni))))
plt.savefig('./Circular/pcolorAlpha.png')
plt.close()

for i in range(len(set(k))):
	kValue = np.array(sorted(list(set(k))))[i]
	plt.plot(Ni[k == kValue], alpha[:, i])
	plt.xlabel('Prey Neighbours', fontsize = 18)
	plt.ylabel(r'$\alpha$', fontsize = 18)
	plt.title(r'$k = '+str(kValue)+'$', fontsize = 20)
	plt.savefig('./Circular/pred'+str(kValue).zfill(3)+'/alpha.png')
