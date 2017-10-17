import numpy as np
import math
import sys
import h5py
import sys
import os
from glob import glob
from scipy.stats import itemfreq
if os.getcwd().rfind('share') > 0:
	import matplotlib as mpl
	mpl.use('Agg')
import matplotlib.pyplot as plt

from matplotlib.colors import LinearSegmentedColormap
from sklearn.neighbors import KDTree
from scipy.spatial import Voronoi
from scipy.spatial.distance import cdist

sys.dont_write_bytecode = True
sys.path.insert(0,os.getcwd())

from defaultParams import *
from fixedParams import *
from params import *

from CSM import init
from CSM import loadData
from CSM import saveData

files = glob('data*h5')

for dataFile in files:
    data = init()
    itime = loadData(data, dataFile)
    saveData(data, int(dataFile[dataFile.rfind('-')+1:dataFile.rfind('.')]), itime)
