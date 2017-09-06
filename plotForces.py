import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import rc
import sys
import h5py
import os
from glob import glob

f, axarr = plt.subplots(4, sharex=True)
size = f.get_size_inches()
f.set_figheight(size[1]*2.)
dataFile = glob('*/*.h5')
save = []

pwd = os.getcwd()
if pwd.rfind('Group') > 0:
    lab = 'Group size = '
elif pwd.rfind('Length') > 0:
    lab = 'Length = '
elif pwd.rfind('Agressiveness') > 0:
    lab = '$p = $'
elif pwd.rfind('noise') > 0:
    lab = '$\eta = $'
if pwd.rfind('predOff') > 0:
    lab = 'Group size = '
    cols = ['blue', 'green', 'red', 'cyan']
else:
    lab = ''

i = -1
for dFile in dataFile:
    i += 1
    if pwd.rfind('Group') > 0:
        value = str(int(dFile[dFile.rfind('p')+1:dFile.rfind('/')].replace('-', '.')))
    elif pwd.rfind('Length') > 0:
        value = dFile[dFile.rfind('th')+2:dFile.rfind('/')].replace('-', '.')
    elif pwd.rfind('Agressiveness') > 0:
        value = dFile[dFile.rfind('p')+1:dFile.rfind('/')].replace('-', '.')
    elif pwd.rfind('noise') > 0:
        value = dFile[dFile.rfind('eta')+3:dFile.rfind('/')].replace('-', '.')
    elif pwd.rfind('predOff') > 0:
        value = str(int(dFile[dFile.rfind('p')+1:dFile.rfind('/')].replace('-', '.')))
        print value
    else:
        value = dFile[:dFile.rfind('/')]

    data = dict()
    h5f = h5py.File(dFile,'r')
    itime = np.copy(h5f['itime'])[0]
    data['forceWalls'] = np.copy(h5f['forceWalls'])
    data['t'] = np.copy(h5f['t'])
    data['sheepAccFlight'] = np.copy(h5f['sheepAccFlight'])
    data['sheepAccFlocking'] = np.copy(h5f['sheepAccFlocking'])
    data['dist_ij'] = np.copy(h5f['dist_ij-mean'])
#    t = ['t_' + dFile[dFile.rfind('th')+2:dFile.rfind('/')].replace('-', '.')]
#    alive = ['alive_' + dFile[dFile.rfind('th')+2:dFile.rfind('/')].replace('-', '.')]
    axarr[0].plot(data['t'][:itime], data['dist_ij'][:itime], lw = 2, label = lab + value, color = cols[i])
    axarr[1].plot(data['t'][:itime], ((data['sheepAccFlocking'][:itime]**2).sum(axis = 2)).mean(axis = 1), lw = 2, color = cols[i])
    axarr[2].plot(data['t'][:itime], ((data['forceWalls'][:itime]**2).sum(axis = 2)).mean(axis = 1), lw = 2, color = cols[i])
    axarr[3].plot(data['t'][:itime], int(value)*((data['sheepAccFlocking'][:itime]**2).sum(axis = 2)).mean(axis = 1), lw = 2, color = cols[i])
#    t = t + map(str, data['t'][:itime])
#    alive = alive + map(str, data['alive'][:itime].sum(axis = 1))
#    save = save + [t] + [alive]

axarr[0].set_ybound(0.6, 1.4)
axarr[1].set_ybound(0,2)
axarr[2].set_ybound(0,0.5)
axarr[3].set_ybound(0,30)

axarr[0].set_ylabel('Separation', fontsize = 18)
axarr[1].set_ylabel('Flocking', fontsize = 18)
axarr[2].set_ylabel('Walls', fontsize = 18)
axarr[3].set_ylabel('Flocking*GroupSize', fontsize = 18)

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

axarr[0].set_xlabel('Time', fontsize = 18)
axarr[1].set_xlabel('Time', fontsize = 18)
axarr[2].set_xlabel('Time', fontsize = 18)
axarr[3].set_xlabel('Time', fontsize = 18)

lgd = axarr[0].legend(loc = 'upper right', bbox_to_anchor=(1.4, 1.))
for label in axarr[0].get_xticklabels() + axarr[0].get_yticklabels() + axarr[1].get_xticklabels() + axarr[1].get_yticklabels() + axarr[2].get_xticklabels() + axarr[2].get_yticklabels():
    label.set_fontsize(16)
plt.savefig('./forces.png', bbox_extra_artists=(lgd,axarr[0],axarr[1],axarr[2],), bbox_inches='tight')
#maxLen = 0
#for i in range(np.shape(save)[0]):
#    if maxLen < len(save[i]):
#        maxLen = len(save[i])

#for i in range(np.shape(save)[0]):
#    while len(save[i]) < maxLen:
#        save[i].append('')

#np.savetxt('./predation.csv', np.transpose(np.array(save)), fmt = '%s', delimiter = ',')
