import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import rc
import sys
import h5py
import os
from glob import glob


dataFile = glob('*/*.h5')
save = []

pwd = os.getcwd()
if pwd.rfind('Group') > 0:
    lab = 'Group size = '
elif pwd.rfind('Length') > 0:
    lab = 'Length = '
else:
    lab = ''

for dFile in dataFile:
    if pwd.rfind('Group') > 0:
        value = str(int(dFile[dFile.rfind('p')+1:dFile.rfind('/')].replace('-', '.')))
    elif pwd.rfind('Length') > 0:
        value = dFile[dFile.rfind('th')+2:dFile.rfind('/')].replace('-', '.')
    else:
        value = dFile[:dFile.rfind('/')]

    data = dict()
    h5f = h5py.File(dFile,'r')
    itime = np.copy(h5f['itime'])[0]
    data['alive'] = np.copy(h5f['alive'])
    data['t'] = np.copy(h5f['t'])
    t = ['t_' + dFile[dFile.rfind('th')+2:dFile.rfind('/')].replace('-', '.')]
    alive = ['alive_' + dFile[dFile.rfind('th')+2:dFile.rfind('/')].replace('-', '.')]
    plt.plot(data['t'][:itime], data['alive'][:itime].sum(axis = 1), lw = 2, label = lab + value)
    t = t + map(str, data['t'][:itime])
    alive = alive + map(str, data['alive'][:itime].sum(axis = 1))
    save = save + [t] + [alive]

plt.ylim(0,405)
plt.axhline(400, color = 'red', ls = '--')
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
plt.xlabel('Time', fontsize = 18)
plt.ylabel('$N(t)$', fontsize = 18)
ax = plt.subplot(111)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(16)
plt.legend(loc = 'upper right')
plt.savefig('./predation.png')
maxLen = 0
for i in range(np.shape(save)[0]):
    if maxLen < len(save[i]):
        maxLen = len(save[i])

for i in range(np.shape(save)[0]):
    while len(save[i]) < maxLen:
        save[i].append('')

np.savetxt('./predation.csv', np.transpose(np.array(save)), fmt = '%s', delimiter = ',')
