import numpy as np
import math
import sys
import h5py
from glob import glob
import os
if os.getcwd().rfind('share'):
	import matplotlib as mpl
	mpl.use('Agg')
import matplotlib.pyplot as plt
sys.dont_write_bytecode = True
sys.path.insert(0,os.getcwd())

from params import *
from defaultParams import *
from fixedParams import *

def initPlot(data, savePlotPng):
	if predOff == False:
		dogTheta = np.arctan2(data['dogVel'][-1,1], data['dogVel'][-1,0])
	sheepTheta = np.arctan2(data['sheepVel'][-1,:,1], data['sheepVel'][-1,:,0])

	from matplotlib import colors
	cmap = colors.ListedColormap(['black', 'blue'])
	bounds=[0,0.5,1.0]
	norm = colors.BoundaryNorm(bounds, cmap.N)
	if predOff == False:
		dogQuiver = plt.quiver(data['dog'][0, 0], data['dog'][0, 1], np.cos(dogTheta), np.sin(dogTheta), scale = 30, color = 'red')
	sheepQuiver = plt.quiver(data['sheep'][0,:,0], data['sheep'][0,:,1], np.cos(sheepTheta), np.sin(sheepTheta), scale = 30, cmap = cmap, norm = norm)

	wallLeft = -wallSize
	wallRight = wallSize
	wallBottom = -wallSize
	wallTop = wallSize
	plt.axis([wallLeft-0.25,wallRight+0.25,wallBottom-0.25,wallTop+0.25])
	plt.axes().set_aspect('equal')
	plotid = 0
	if wallType == 'Square':
		plt.axhline(wallTop, color = 'r', lw = 5)
		plt.axhline(wallBottom, color = 'r', lw = 5)
		plt.axvline(wallLeft, color = 'r', lw = 5)
		plt.axvline(wallRight, color = 'r', lw = 5)
	elif wallType == 'Circular':
		circle = plt.Circle((0,0), data['walls'][0], color = 'r', lw = 5, fill = False)
		fig = plt.gcf()
		ax = fig.gca()
		ax.add_artist(circle)

	if savePlotPng == 'On':
		plt.savefig('frames/'+str(0).zfill(7)+'.png')
	else:
		plt.pause(0.005)

	if predOff == False:
		return dogQuiver, sheepQuiver
	else:
		return sheepQuiver

def plotDataPositions(data, tstep, dQuiv, sQuiv, savePlotPng):
	if predOff == False:
		dogTheta = np.arctan2(data['dogVel'][tstep-1,1], data['dogVel'][tstep-1,0])
	if constantNP == False:
		sheepTheta = np.arctan2(data['sheepVel'][tstep-1,:,1][data['alive'][tstep]], data['sheepVel'][tstep-1,:,0][data['alive'][tstep]])
	else:
		sheepTheta = np.arctan2(data['sheepVel'][tstep-1,:,1], data['sheepVel'][tstep-1,:,0])

	colorQuiver = np.zeros(NP)
	if predOff == False:
		if showSheepDogCanSee == 'On':
			colorQuiver[data['dog_index_neighbours']] = 1.
			if showDogInfluence == 'On' or segmentColours == 'On' or constantNPShow == 'On':
				sys.exit('You cannot show multiple colourings')
		elif showDogInfluence == 'On':
			colorQuiver[(sheepTheta > dogTheta - np.pi/16 ) * (sheepTheta < dogTheta + np.pi/16)] = 1
			if showSheepDogCanSee == 'On' or segmentColours == 'On' or constantNPShow == 'On':
				sys.exit('You cannot show multiple colourings')
		elif segmentColours == 'On':
			colorQuiver = data['interactingSheep'][tstep-1]
			if showDogInfluence == 'On' or showSheepDogCanSee == 'On' or constantNPShow == 'On':
				sys.exit('You cannot show multiple colourings')
		elif constantNPShow == 'On':
			colorQuiver[data['alive'][tstep-1] == False] = 1.
			if showDogInfluence == 'On' or showSheepDogCanSee == 'On' or segmentColours == 'On':
				sys.exit('You cannot show multiple colourings')

		dQuiv.set_offsets(np.transpose([data['dog'][tstep, 0], data['dog'][tstep, 1]]))
		dQuiv.set_UVC(np.cos(dogTheta),np.sin(dogTheta))
	if constantNP == False:
		sQuiv.set_offsets(np.transpose([data['sheep'][tstep,:, 0][data['alive'][tstep]], data['sheep'][tstep,:, 1][data['alive'][tstep]]]))
		sQuiv.set_UVC(np.cos(sheepTheta), np.sin(sheepTheta), colorQuiver)
	else:
		sQuiv.set_offsets(np.transpose([data['sheep'][tstep,:, 0], data['sheep'][tstep,:, 1]]))
		sQuiv.set_UVC(np.cos(sheepTheta), np.sin(sheepTheta), colorQuiver)

	if savePlotPng == 'On':
		plt.savefig('frames/'+str(int(np.floor(data['t'][tstep]/plotPeriod))).zfill(7)+'.png')
	else:
		plt.pause(0.005)

files = sorted(glob('*-1.h5'))

for dataName in files:
	data = dict()
	h5f = h5py.File(dataName,'r')
	itime = np.copy(h5f['itime'])[0]

	if dataName == files[0]:
		data['t'] = np.copy(h5f['t'])[0:itime]
		data['sheep'] = np.copy(h5f['sheep'])[0:itime]
		data['sheepVel'] = np.copy(h5f['sheepVel'])[0:itime]
		data['dog'] = np.copy(h5f['dog'])[0:itime]
		data['dogVel'] = np.copy(h5f['dogVel'])[0:itime]
		if segments == 'On':
			data['interactingSheep'] = data['interactingSheep'][0:itime]
	else:
		data['t'] = np.copy(h5f['t'])[4:itime]
		data['sheep'] = np.copy(h5f['sheep'])[4:itime]
		data['sheepVel'] = np.copy(h5f['sheepVel'])[4:itime]
		data['dog'] = np.copy(h5f['dog'])[4:itime]
		data['dogVel'] = np.copy(h5f['dogVel'])[4:itime]
		if segments == 'On':
			data['interactingSheep'] = data['interactingSheep'][4:itime]

	tstep = 0
	lastPlot = data['t'][tstep]
	if dataName == files[0]:
		if predOff == False:
			dogQuiver, sheepQuiver = initPlot(data, savePlotPng)
		else:
			sheepQuiver = initPlot(data, savePlotPng)
	for t in data['t']:
		if t-lastPlot > plotPeriod:
			print 'read', t
			if predOff == False:
				plotDataPositions(data, tstep, dogQuiver, sheepQuiver, savePlotPng)
			else:
				plotDataPositions(data, tstep, 'off', sheepQuiver, savePlotPng)
			lastPlot = data['t'][tstep]
		tstep += 1
