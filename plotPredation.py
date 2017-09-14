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


dataFile = glob('*/*-1.h5')
save = []

def func(x, a, b):
    return a*np.exp(-b*x)

pwd = os.getcwd()
if pwd.rfind('Group') > 0:
    lab = 'Group size = '
elif pwd.rfind('Length') > 0:
    lab = 'Length = '
elif pwd.rfind('Agressiveness') > 0:
    lab = '$p = $'
elif pwd.rfind('noise') > 0:
    lab = '$\eta = $'
else:
    lab = ''

a,b,G = [[],[],[]]
for dFile in dataFile:
    if pwd.rfind('Group') > 0:
        value = str(int(dFile[dFile.rfind('p')+1:dFile.rfind('/')].replace('-', '.')))
    elif pwd.rfind('Length') > 0:
        value = dFile[dFile.rfind('th')+2:dFile.rfind('/')].replace('-', '.')
    elif pwd.rfind('Agressiveness') > 0:
        value = dFile[dFile.rfind('p')+1:dFile.rfind('/')].replace('-', '.')
    elif pwd.rfind('noise') > 0:
        value = dFile[dFile.rfind('eta')+3:dFile.rfind('/')].replace('-', '.')
    else:
        value = dFile[:dFile.rfind('/')]

    data = dict()
    h5f = h5py.File(dFile,'r')
    itime = np.copy(h5f['itime'])[0]
    data['alive'] = np.copy(h5f['alive'])
    data['t'] = np.copy(h5f['t'])

    if dFile == dataFile[0]:
        NP = np.shape(data['alive'])[1]

    t = ['t_' + dFile[dFile.rfind('th')+2:dFile.rfind('/')].replace('-', '.')]
    alive = ['alive_' + dFile[dFile.rfind('th')+2:dFile.rfind('/')].replace('-', '.')]
    plt.semilogy(data['t'][:itime], data['alive'][:itime].sum(axis = 1), lw = 2, label = lab + value)
    t = t + map(str, data['t'][:itime])
    alive = alive + map(str, data['alive'][:itime].sum(axis = 1))
    save = save + [t] + [alive]
    popt, pcov = curve_fit(func, data['t'][:itime], data['alive'][:itime].sum(axis = 1))
    a = a + [popt[0]]
    b = b + [popt[1]]
    G = G + [int(value)]
    print lab + value, a[-1], b[-1]


plt.axhline(np.log(NP), color = 'red', ls = '--')
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
plt.xlabel('Time', fontsize = 18)
plt.ylabel('$\log(N(t))$', fontsize = 18)
ax = plt.subplot(111)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(16)
plt.legend(loc = 'lower left')
plt.savefig('./logYPredation.png')

plt.close()


for dFile in dataFile:
    if pwd.rfind('Group') > 0:
        value = str(int(dFile[dFile.rfind('p')+1:dFile.rfind('/')].replace('-', '.')))
    elif pwd.rfind('Length') > 0:
        value = dFile[dFile.rfind('th')+2:dFile.rfind('/')].replace('-', '.')
    elif pwd.rfind('Agressiveness') > 0:
        value = dFile[dFile.rfind('p')+1:dFile.rfind('/')].replace('-', '.')
    elif pwd.rfind('noise') > 0:
        value = dFile[dFile.rfind('eta')+3:dFile.rfind('/')].replace('-', '.')
    else:
        value = dFile[:dFile.rfind('/')]

    data = dict()
    h5f = h5py.File(dFile,'r')
    itime = np.copy(h5f['itime'])[0]
    data['alive'] = np.copy(h5f['alive'])
    data['t'] = np.copy(h5f['t'])
    plt.plot(data['t'][:itime], data['alive'][:itime].sum(axis = 1), lw = 2, label = lab + value)

plt.ylim(0,NP + 5)
plt.xlim(0,100)
plt.axhline(NP, color = 'red', ls = '--')
if topsy == False:
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)
plt.xlabel('Time', fontsize = 18)
plt.ylabel('$N(t)$')
plt.legend(loc = 'upper right', fontsize = 16)
plt.savefig('./predation.png')

plt.close()

plt.scatter(G, a, label = 'a', color = 'b')
plt.scatter(G, b, label = 'b', color = 'g')
plt.ylabel('Value', fontsize = 18)
plt.xlabel('Group Size', fontsize = 18)
plt.legend(loc = 'upper right')
plt.savefig('./fitValues.png')

plt.close()

plt.scatter(G, a, label = 'a')
plt.ylabel('Value', fontsize = 18)
plt.xlabel('Group Size', fontsize = 18)
plt.savefig('./fitValuesA.png')

plt.close()

plt.scatter(G, b, label = 'b')
plt.ylabel('Value', fontsize = 18)
plt.xlabel('Group Size', fontsize = 18)
plt.savefig('./fitValuesB.png')

plt.close()

for i in range(len(dataFile)):
    dFile = dataFile[i]
    data = dict()
    h5f = h5py.File(dFile,'r')
    itime = np.copy(h5f['itime'])[0]
    data['alive'] = np.copy(h5f['alive'])
    data['t'] = np.copy(h5f['t'])
    plt.semilogy(np.linspace(0, 100, 1000), func(np.linspace(0, 100, 1000), a[i], b[i]), label = 'Fit')
    plt.semilogy(data['t'][:itime][0::100], data['alive'][:itime].sum(axis = 1)[0::100], label = 'Data', color = 'k', marker = 'o')
    plt.legend(loc = 'lower left')
    plt.xlabel('Time', fontsize = 18)
    plt.ylabel('$N(t)$', fontsize = 18)
    ax = plt.subplot(111)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(16)
    plt.savefig('./fit' + dFile[:dFile.rfind('/')] + '.png')
    plt.close()

maxLen = 0
for i in range(np.shape(save)[0]):
    if maxLen < len(save[i]):
        maxLen = len(save[i])

for i in range(np.shape(save)[0]):
    while len(save[i]) < maxLen:
        save[i].append('')

np.savetxt('./predation.csv', np.transpose(np.array(save)), fmt = '%s', delimiter = ',')
