import numpy as np
import matplotlib.pyplot as plt
from glob import glob

volumeFiles = glob('*/vol*')

area = []
for file in volumeFiles:
    openFile = open(file, 'r')
    area.append(openFile.read())
    openFile.close()

n = []
for file in volumeFiles:
    n.append(int(file[5:-7]))

plt.plot(n, area, 'o')
plt.ylabel('Area')
plt.xlabel('Size of interacting group')
plt.xlim(-0.1, 405)
plt.show()
#plt.safefig('volumeVSgroupSize.png')
