import numpy as np
import h5py
import matplotlib.pyplot as plt
from glob import glob
from scipy.optimize import curve_fit

def func(x, a, b):
    return a*x + b#a*x**(-b)

dataFiles = glob('*/*h5')

m1 = []
m2 = []
Ni = []

for dFile in dataFiles:
    data = dict()
    h5f = h5py.File(dFile,'r')
    itime = np.copy(h5f['itime'])[0]
    data['sheep'] = np.copy(h5f['sheep'])

    Ni = Ni + [int(dFile[dFile.rfind('p')+1:dFile.rfind('/')].replace('-', '.'))]

    m2 = m2 + [np.sort(np.sqrt((data['sheep'][itime]**2).sum(axis = 1)))[-50:].mean()]
    m1 = m1 + [np.sort(np.sqrt((data['sheep'][itime]**2).sum(axis = 1)))[10:].mean()]

    #file = open(dFile[:dFile.rfind('/')]+'/output/volume', 'r')
    #volume = float(file.read())
    #m2 = m2 + [np.sqrt(volume/np.pi)]


poptOutter, pcovOutter = curve_fit(func, Ni, m2)
print 'Outter:', poptOutter

poptInner, pcovInner = curve_fit(func, Ni, m1)
print 'Inner:',poptInner

plt.plot(np.linspace(0.01, 500, 1000), func(np.linspace(0.01, 500, 1000), poptOutter[0], poptOutter[1]), label = 'Fit', color = 'k')
plt.plot(np.linspace(0.01, 500, 1000), func(np.linspace(0.01, 500, 1000), poptInner[0], poptInner[1]), color = 'k')
plt.scatter(Ni, m1, label = 'Inner Radius')
plt.scatter(Ni, m2, label = 'Outer Radius')
plt.legend(loc = 'upper right')
plt.xlabel('Group size, $Ni$')
plt.ylabel('Radius of steady state, $R$')
plt.ylim(0, 2)
plt.savefig('radiusNi.png')
