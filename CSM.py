import numpy as np
import math
import sys
import h5py
import sys
import os
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

from CSMfunctions.noiseFunc import *
from CSMfunctions.flockingFunc import *
from CSMfunctions.flightFunc import *
from CSMfunctions.wallFunc import *

def makeSquareWalls(wallSize):
	#Walls defined by four lines: Ax+By+C = 0
	w = dict()
	w['eqn'] = [[0,1,-wallSize],[0,1,wallSize],[1,0,wallSize],[1,0,-wallSize]]
	w['n'] = [np.array([0,-1]),np.array([0,1]),np.array([1,0]),np.array([-1,0])]
	return w

def init():
	#Sets up main data structure
	#Returns: dictionary of simulation data
	data = dict()
	cachedTimesteps = 6*int(round(snapshotPeriod/dt))
	data['dog'] = np.zeros((cachedTimesteps, 2))
	data['dogVel'] = np.zeros((cachedTimesteps,2))
	data['sheep'] = np.zeros((cachedTimesteps,NP, 2))
	data['sheepVel'] = np.zeros((cachedTimesteps,NP,2))

	data['t'] = np.zeros(cachedTimesteps)
	data['dog_index_neighbours'] = []

	if timeStepMethod != 'Euler' and timeStepMethod != 'Adaptive':
		sys.exit("Error: updateMethod not recognized!")

	if wallType == 'Square' :
		data['walls'] = makeSquareWalls(wallSize)
		data['forceWalls'] = np.zeros((cachedTimesteps,NP,2))
		data['forceWallsDog'] = np.zeros((cachedTimesteps,2))
	elif wallType == 'Circular':
		data['walls'] = [wallSize]
		data['forceWalls'] = np.zeros((cachedTimesteps,NP,2))
		data['forceWallsDog'] = np.zeros((cachedTimesteps,2))
	else:
		data['walls'] = []

	if sheepSheepInteration == 'On':
		data['forceSheep'] = np.zeros((cachedTimesteps, NP, 2))

	data['alive'] = np.zeros((cachedTimesteps, NP)) == 0

	if updateMethod == 'Velocity':
		data['bterm'] = np.zeros((cachedTimesteps,NP,2))
		data['aterm'] = np.zeros((cachedTimesteps,NP,2))
	elif updateMethod == 'Acceleration':
		data['sheepAccFlight'] = np.zeros((cachedTimesteps,NP,2))
		data['sheepAcc'] = np.zeros((cachedTimesteps,NP,2))
		data['dogAcc'] = np.zeros((cachedTimesteps,2))
	else:
		sys.exit("Error: updateMethod not recognized!")
	if segments == 'On':
		data['interactingSheep'] = np.zeros((cachedTimesteps, NP))
	if constantNP == True:
		data['deathCount'] = np.zeros(cachedTimesteps)
	return data

def initCond(data):
	data['dog'][0] = np.array(dog_init)
	data['dogVel'][0] = np.array(dog_vel_init)
	if sheep_init == 'Grid':
		data['sheep'][0] = 2*np.array([[x,y] for x in np.arange(0.,np.ceil(np.sqrt(NP)),1.) for y in np.arange(0.,np.ceil(np.sqrt(NP)),1.)])[0:NP]
		data['sheepVel'][0] = np.zeros((NP, 2))
	elif sheep_init == 'Random':
		print 'Creating initial locations ...'
		data['sheep'][0,:,0] = np.random.rand(NP)*sheep_area - sheep_area/2.
		data['sheep'][0,:,1] = np.random.rand(NP)*sheep_area - sheep_area/2.
		sheepTheta = np.random.rand(NP)*2*math.pi

		data['sheepVel'][0,:,0] = sheep_vel_init*np.cos(sheepTheta)
		data['sheepVel'][0,:,1] = sheep_vel_init*np.sin(sheepTheta)
		data['sheepVel'][-1,:,0] = sheep_vel_init*np.cos(sheepTheta)
		data['sheepVel'][-1,:,1] = sheep_vel_init*np.sin(sheepTheta)
	elif sheep_init == 'Anular':
		posPrey = np.linspace(0, 2*np.pi, NP)
		data['sheep'][0,:,0] = np.cos(posPrey)
		data['sheep'][0,:,1] = np.sin(posPrey)


###############
#Velocity
###############
def doVelocityStep(data):
	#Velocity of sheep
	distMatrixNotMe = np.array(map(lambda j:data['sheep'][itime,j] - np.delete(data['sheep'][itime], j, 0),range(NP)))
	normStuff = np.transpose(np.tile(np.transpose(np.linalg.norm(distMatrixNotMe, axis = 2)**(-2)),(2,1,1)))
	preyMinusPred = data['sheep'][itime] - data['dog'][itime]
	normstuff2 = np.transpose(np.tile(np.linalg.norm(preyMinusPred,axis=1),(2,1)))
	data['bterm'][itime] = b*preyMinusPred*normstuff2**(-2)
	data['sheepVel'][itime] = (distMatrixNotMe*(normStuff - a)).sum(axis = 1)/NP + data['bterm'][itime]

	#Velocity of dog
	distDogSheep = data['sheep'][itime] - data['dog'][itime]
	frac = distDogSheep*np.transpose(np.tile((np.linalg.norm(distDogSheep,axis=1))**(-p),(2,1)))
	data['dogVel'][itime] = (c/NP)*frac.sum(axis=0)

	if timeStepMethod == 'Euler':
		data['dog'][itime + 1] = data['dog'][itime] + data['dogVel'][itime]*dt
		data['sheep'][itime + 1] = data['sheep'][itime] + data['sheepVel'][itime]*dt
		data['t'][itime+1] = data['t'][itime] + dt
	elif timeStepMethod == 'Adaptive':
		if itime == 0:
			data['sheep'][itime + 1] = data['sheep'][itime] + data['sheepVel'][itime]*dt
			data['dog'][itime + 1] = data['dog'][itime] + data['dogVel'][itime]*dt
			data['t'][itime+1] = data['t'][itime] + dt
		elif itime == 1:
			data['sheep'][itime + 1] = data['sheep'][itime] + 1.5*dt*data['sheepVel'][itime] - 0.5*dt*data['sheepVel'][itime - 1]
			data['dog'][itime +1] = data['dog'][itime] + data['dogVel'][itime]*dt
			data['t'][itime+1] = data['t'][itime] + dt
		elif itime == 2:
			data['sheep'][itime + 1] = data['sheep'][itime] + dt*((23./12.)*data['sheepVel'][itime] - (4./3.)*data['sheepVel'][itime - 1] + (5./12.)*data['sheepVel'][itime - 2])
			data['dog'][itime + 1] = data['dog'][itime] + data['dogVel'][itime]*dt
			data['t'][itime+1] = data['t'][itime] + dt
		elif itime == 3:
			sheepTemp = data['sheep'][itime] + (dt/24.)*(55.*data['sheepVel'][itime] - 59.*data['sheepVel'][itime - 1] + 37.*data['sheepVel'][itime - 2] - 9.*data['sheepVel'][itime - 3])
			sheepVelTemp = (sheepTemp - data['sheep'][itime])/dt
			data['sheep'][itime + 1] = data['sheep'][itime] + (dt/24.)*(9.*sheepVelTemp + 19.*data['sheepVel'][itime] - 5.*data['sheepVel'][itime - 1] + data['sheepVel'][itime - 2])
			data['dog'][itime + 1] = data['dog'][itime] + data['dogVel'][itime]*dt
			q_f = np.min(((270./19.)*((epsilon*dt)/(np.abs(data['sheep'][itime + 1] - sheepTemp))))**(0.25))
			q = np.floor(100.*q_f)/100.
			data['t'][itime+1] = data['t'][itime] + dt
		else:
			sheepTemp = data['sheep'][itime] + q*(dt/24.)*(55.*data['sheepVel'][itime] - 59.*data['sheepVel'][itime - 1] + 37.*data['sheepVel'][itime - 2] - 9.*data['sheepVel'][itime - 3])
			sheepVelTemp = (sheepTemp - data['sheep'][itime])/(dt*q)
			data['sheep'][itime + 1] = data['sheep'][itime] + q*(dt/24.)*(9.*sheepVelTemp + 19.*data['sheepVel'][itime] - 5.*data['sheepVel'][itime - 1] + data['sheepVel'][itime - 2])
			data['dog'][itime + 1] = data['dog'][itime] + data['dogVel'][itime]*dt*q
			q_f = np.min(((270./19.)*((epsilon*dt)/(np.abs(data['sheep'][itime + 1] - sheepTemp))))**(0.25))
			data['t'][itime+1] = data['t'][itime] + dt*q
			q = np.floor(100.*q_f)/100.

	else:
		sys.exit('Invalid time step method:'+timeStepMethod+' in doVelocityStep()')

###############
#Acceleration
###############
def doAccelerationStep(data, q=0):
	#Prints the location of the prey furthest out
	if emergencyCheck == 'On':
		print 'furthest up', data['sheep'][itime,np.where(np.max(data['sheep'][itime][:,1]) == data['sheep'][itime][:,1]),:]
		print 'furthest right', data['sheep'][itime,np.where(np.max(data['sheep'][itime][:,0]) == data['sheep'][itime][:,0]),:]
		print 'furthest down',data['sheep'][itime,np.where(np.min(data['sheep'][itime][:,1]) == data['sheep'][itime][:,1]),:]
		print 'furthest left', data['sheep'][itime,np.where(np.min(data['sheep'][itime][:,0]) == data['sheep'][itime][:,0]),:]


	dist_id = cdist(data['sheep'][itime], data['dog'][itime].reshape(1,2)).reshape(NP)
	if predation == 'On':
		if np.min(dist_id[data['alive'][itime]]) < predation_length_scale:
			prey_too_close = np.where(dist_id == np.min(dist_id[data['alive'][itime]]))[0][0]

			if constantNP == False:
				print 'predation of prey', prey_too_close, ': ', sum(data['alive'][itime])
				data['alive'][itime:, prey_too_close] = False
			else:
				print 'Moving prey: ', prey_too_close
				min_dist_ij = 0.

				while min_dist_ij < 0.001:
					data['sheep'][itime, prey_too_close] = np.random.rand(2)*wallSize - wallSize/2
					dist_ij = cdist(data['sheep'][itime, prey_too_close].reshape(1,2), data['sheep'][itime])
					min_dist_ij = np.min(dist_ij[np.nonzero(dist_ij)])
				data['sheepVel'][itime, prey_too_close] = 0.0
				data['sheepAcc'][itime, prey_too_close] = 0.0
				data['deathCount'][itime:] = data['deathCount'][itime - 1] + 1


	data = includeNoise(data, itime, noise)

	#Acceleration of sheep
	#Flocking
	if flocking == 'PredatorPrey':
		sheepAccFlocking = doPredatorPrey(data, itime, NP, f, a, groupSize, sheepSize, predation, gaussian)

	elif flocking == 'PPPoisson':
		sheepAccFlocking = doPPPoisson(data, itime, NP, f, a, lam, groupSize, predation)

	elif flocking == 'Vicsek':
		sheepAccFlocking = doVicsek(data, itime, NP, r)

	elif flocking == 'Topo':
		sheepAccFlocking = doTopo(data, itime, NP, n)

	else:
		sys.exit('Invalid flocking mechanisim: ' + flocking + ' in doAccelerationStep')

	#Flight
	data = doFlight(data, itime, NP, b, allSheepSeeDog)

	data['sheepAcc'][itime][data['alive'][itime]] = (data['sheepAccFlight'][itime][data['alive'][itime]] + sheepAccFlocking - data['sheepVel'][itime][data['alive'][itime]])/sheepMass

	if cap == 'On':
		indices = np.where(np.sqrt((data['sheepAcc'][itime]**2).sum(axis = 1)) > accCap)
		data['sheepAcc'][itime][indices] = accCap*data['sheepAcc'][itime][indices]/np.transpose(np.tile(np.sqrt((data['sheepAcc'][itime][indices]**2).sum(axis = 1)), (2, 1)))

	data = doWalls(data, itime, A, B, C, NP, wallType)


	data['sheepAcc'][itime][data['alive'][itime]] += data['forceWalls'][itime][data['alive'][itime]]

#################################################################
	#Acceleration of dog
	if segments == 'On':
		angles = [np.arctan2(i[1], i[0]) for i in np.tile(data['dog'][itime], (NP, 1)) - data['sheep'][itime][data['alive'][itime]]]
		freq,locs = np.histogram(angles, bins = noSegments, density = False, normed = False)
		maxLoc = np.argmax(freq)
		rangeLoc = locs[maxLoc: maxLoc+2]
		sheep = data['sheep'][itime][data['alive'][itime]][(np.array(angles) > rangeLoc[0]) &  (np.array(angles) < rangeLoc[1])]
		numberInteractingSheep = np.shape(sheep)[0]
		data['interactingSheep'][itime] = ((np.array(angles) > rangeLoc[0]) &  (np.array(angles) < rangeLoc[1])).tolist()
	elif dogNeareastNeigh == 'On':
		tree = KDTree(data['sheep'][itime][data['alive'][itime]])
		idx = tree.query(data['dog'][itime].reshape(1, -1), np.min([k, sum(data['alive'][itime])]))[1][0]
		sheep = data['sheep'][itime][data['alive'][itime]][idx]
		numberInteractingSheep = k
	else:
		sheep = data['sheep'][itime][data['alive'][itime]]
		numberInteractingSheep = sum(data['alive'][itime])

	preyMinusPred = sheep - data['dog'][itime]

	normStuff3 = np.linalg.norm(preyMinusPred,axis=1).reshape(numberInteractingSheep,1)
	Hfunc = np.zeros((NP, 2))
	Hfunc[0:numberInteractingSheep] = 1./(normStuff3**p + 0.1)
	Hfunc[numberInteractingSheep:] = np.nan
	data['dogAcc'][itime] = (-data['dogVel'][itime] + c/numberInteractingSheep*(preyMinusPred*(normStuff3**(-1))*Hfunc[:numberInteractingSheep]).sum(axis=0))/dogMass


	if normDog == 'On':
		data['dogAcc'][itime] = data['dogAcc'][itime]/np.transpose(np.tile(np.sqrt((data['dogAcc'][itime]**2).sum()) ,(2,1)))

	if wallType == 'Square':
		for wall in range(len(data['walls']['eqn'])):
			w = data['walls']['eqn'][wall]

			data['forceWallsDog'][itime] += (A*np.exp((np.dot(-data['walls']['n'][wall], data['dog'][itime].T)-np.abs(w[2])+C)/B)*data['walls']['n'][wall])
	elif wallType == 'Circular':
		distance = np.linalg.norm(data['dog'][itime])
		if distance == 0.0:
			distance = 0.01
		unit = data['dog'][itime]/distance

		data['forceWallsDog'][itime] = -A*unit*np.exp((distance-data['walls'][0]+C)/B )


	data['dogAcc'][itime] = data['dogAcc'][itime] + data['forceWallsDog'][itime]




##############################################################################
	if timeStepMethod == 'Euler':
		#Euler time step
		data['dogVel'][itime + 1] = data['dogVel'][itime] + data['dogAcc'][itime]*dt
		data['sheepVel'][itime + 1] = data['sheepVel'][itime] + data['sheepAcc'][itime]*dt
		data['dog'][itime + 1] = data['dog'][itime] + data['dogVel'][itime]*dt
		data['sheep'][itime + 1] = data['sheep'][itime] + data['sheepVel'][itime]*dt
		data['t'][itime+1] = data['t'][itime] + dt
		return 0
	elif timeStepMethod == 'Adaptive':
		#Adaptive time step
		if itime == 0:
			data['dogVel'][itime + 1] = data['dogVel'][itime] + data['dogAcc'][itime]*dt
			data['sheepVel'][itime + 1] = data['sheepVel'][itime] + data['sheepAcc'][itime]*dt
			data['sheep'][itime + 1] = data['sheep'][itime] + data['sheepVel'][itime]*dt
			data['dog'][itime + 1] = data['dog'][itime] + data['dogVel'][itime]*dt
			data['t'][itime+1] = data['t'][itime] + dt
			return 0
		elif itime == 1:
			data['sheepVel'][itime + 1] = data['sheepVel'][itime] + 1.5*dt*data['sheepAcc'][itime] - 0.5*dt*data['sheepAcc'][itime - 1]
			data['dogVel'][itime + 1] = data['dogVel'][itime] + 1.5*dt*data['dogAcc'][itime] - 0.5*dt*data['dogAcc'][itime - 1]
			data['sheep'][itime + 1] = data['sheep'][itime] + data['sheepVel'][itime]*dt
			data['dog'][itime + 1] = data['dog'][itime] + data['dogVel'][itime]*dt
			data['t'][itime+1] = data['t'][itime] + dt
			return 0
		elif itime == 2:
			data['sheepVel'][itime + 1] = data['sheepVel'][itime] + dt*((23./12.)*data['sheepAcc'][itime] - (4./3.)*data['sheepAcc'][itime - 1] + (5./12.)*data['sheepAcc'][itime - 2])
			data['dogVel'][itime + 1] = data['dogVel'][itime] + dt*((23./12.)*data['dogAcc'][itime] - (4./3.)*data['dogAcc'][itime - 1] + (5./12.)*data['dogAcc'][itime - 2])
			data['sheep'][itime + 1] = data['sheep'][itime] + data['sheepVel'][itime]*dt
			data['dog'][itime + 1] = data['dog'][itime] + data['dogVel'][itime]*dt
			data['t'][itime+1] = data['t'][itime] + dt
			return 0
		elif itime == 3:
			sheepVelTemp = data['sheepVel'][itime] + (dt/24.)*(55.*data['sheepAcc'][itime] - 59.*data['sheepAcc'][itime - 1] + 37.*data['sheepAcc'][itime - 2] - 9.*data['sheepAcc'][itime - 3])
			sheepAccTemp = (sheepVelTemp - data['sheepVel'][itime])/dt
			data['sheepVel'][itime + 1] = data['sheepVel'][itime] + (dt/24.)*(9.*sheepAccTemp + 19.*data['sheepAcc'][itime] - 5.*data['sheepAcc'][itime - 1] + data['sheepAcc'][itime - 2])
			data['sheep'][itime + 1] = data['sheep'][itime] + data['sheepVel'][itime]*dt

			dogVelTemp = data['dogVel'][itime] + (dt/24.)*(55.*data['dogAcc'][itime] - 59.*data['dogAcc'][itime - 1] + 37.*data['dogAcc'][itime - 2] - 9.*data['dogAcc'][itime - 3])
			dogAccTemp = (dogVelTemp - data['dogVel'][itime])/dt
			data['dogVel'][itime + 1] = data['dogVel'][itime] + (dt/24.)*(9.*dogAccTemp + 19.*data['dogAcc'][itime] - 5.*data['dogAcc'][itime - 1] + data['dogAcc'][itime - 2])
			data['dog'][itime + 1] = data['dog'][itime] + data['dogVel'][itime]*dt

			q_f = np.min(((270./19.)*((epsilon*dt)/(np.abs(data['sheepVel'][itime + 1] - sheepVelTemp))))**(0.25))
			q_d = np.min(((270./19.)*((epsilon*dt)/(np.abs(data['dogVel'][itime + 1] - dogVelTemp))))**(0.25))

			min_q = np.min([q_d, q_f])

			q = np.floor(100.*min_q)/100.
			return q
			data['t'][itime+1] = data['t'][itime] + dt
		else:
			sheepVelTemp = data['sheepVel'][itime] + q*(dt/24.)*(55.*data['sheepAcc'][itime] - 59.*data['sheepAcc'][itime - 1] + 37.*data['sheepAcc'][itime - 2] - 9.*data['sheepAcc'][itime - 3])
			sheepAccTemp = (sheepVelTemp - data['sheepVel'][itime])/(dt*q)
			data['sheepVel'][itime + 1] = data['sheepVel'][itime] + q*(dt/24.)*(9.*sheepAccTemp + 19.*data['sheepAcc'][itime] - 5.*data['sheepAcc'][itime - 1] + data['sheepAcc'][itime - 2])
			data['sheep'][itime + 1] = data['sheep'][itime] + data['sheepVel'][itime]*dt*q

			dogVelTemp = data['dogVel'][itime] + q*(dt/24.)*(55.*data['dogAcc'][itime] - 59.*data['dogAcc'][itime - 1] + 37.*data['dogAcc'][itime - 2] - 9.*data['dogAcc'][itime - 3])
			dogAccTemp = (dogVelTemp - data['dogVel'][itime])/(dt*q)
			data['dogVel'][itime + 1] = data['dogVel'][itime] + q*(dt/24.)*(9.*dogAccTemp + 19.*data['dogAcc'][itime] - 5.*data['dogAcc'][itime - 1] + data['dogAcc'][itime - 2])
			data['dog'][itime + 1] = data['dog'][itime] + data['dogVel'][itime]*dt*q

			q_f = np.min(((270./19.)*((epsilon*dt)/(np.abs(data['sheepVel'][itime + 1] - sheepVelTemp))))**(0.25))
			q_d = np.min(((270./19.)*((epsilon*dt)/(np.abs(data['dogVel'][itime + 1] - dogVelTemp))))**(0.25))

			min_q = np.min([q_d, q_f])

			data['t'][itime+1] = data['t'][itime] + dt*q
			q = np.floor(100.*min_q)/100.
			return q
	else:
		sys.exit('Invalid time step method:'+timeStepMethod+' in doAccelerationStep()')
	if dogCentrePrey == 'On':
		data['dog'][itime+1] = data['sheep'][itime+1][data['alive'][itime]].mean(axis = 0)



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

def saveData(data, ens, step):
	h5f = h5py.File('data-%05d-'%int(data['t'][step])+str(ens)+'.h5', 'w')
	h5f.create_dataset('itime', data = [step])
	for k in data.keys():
		if type(data[k]) is dict:
			for k2 in data[k].keys():
				try:
					h5f.create_dataset(k+'-'+k2, data=data[k][k2][0:step+1],compression="gzip")
				except:
					h5f.create_dataset(k+'-'+k2, data=data[k][k2][0:step+1])
		else:
			try:
				h5f.create_dataset(k, data=data[k][0:step+1],compression="gzip")
			except:
				h5f.create_dataset(k, data=data[k][0:step+1])
	h5f.close()

def loadData(data,fn):
	h5f = h5py.File(fn,'r')
	itime = np.copy(h5f['itime'])[0]
	data['dog'][0:itime+1] = np.copy(h5f['dog'])
	data['sheep'][0:itime+1] = np.copy(h5f['sheep'])
	data['sheepVel'][0:itime+1] = np.copy(h5f['sheepVel'])
	data['dogVel'][0:itime+1] = np.copy(h5f['dogVel'])


	data['t'][0:itime+1] = np.copy(h5f['t'])
	data['alive'][0:itime+1] = np.copy(h5f['alive'])

	if timeStepMethod != 'Euler' and timeStepMethod != 'Adaptive':
		sys.exit("Error: updateMethod not recognized!")

	if wallType == 'Square' :
		data['walls']['eqn'] = np.copy(h5f['walls-eqn'])
		data['walls']['n'] = np.copy(h5f['walls-n'])
		data['forceWalls'][0:itime+1] = np.copy(h5f['forceWalls'])
	if wallType == 'Circular':
		data['walls'] = np.copy(h5f['walls']).tolist()
	else:
		data['walls'] = []

	if updateMethod == 'Velocity':
		data['bterm'][0:itime+1] = np.copy(h5f['bterm'])
		data['aterm'][0:itime+1] = np.copy(h5f['aterm'])
	elif updateMethod == 'Acceleration':
		data['sheepAccFlight'][0:itime+1] = np.copy(h5f['sheepAccFlight'])
		data['sheepAcc'][0:itime+1] = np.copy(h5f['sheepAcc'])
		data['dogAcc'][0:itime+1] = np.copy(h5f['dogAcc'])
	else:
		sys.exit("Error: updateMethod not recognized!")

	if segments == 'On':
		data['interactingSheep'] = np.copy(h5f['interactingSheep'])
	return itime

def wipeData(data):
	for k in data.keys():
		if k != 'walls':
			if k == 'dog_index_neighbours':
				data['dog_index_neighbours'] = []
			elif type(data[k]) is dict:
				for k2 in data[k].keys():
					data[k][k2][0:5] = data[k][k2][itime-4:itime+1]
					data[k][k2][5:itime + 1] = data[k][k2][5:itime + 1]*0.
			else:
				data[k][0:5] = data[k][itime-4:itime+1]
				data[k][5:itime + 1] = data[k][5:itime + 1]*0.


###############
#Simulation
###############
def main():
	pass

if __name__ == "__main__":
	plt.close()
	data = init()
	itime = 0
	initCond(data)
	if len(sys.argv) > 1:
		e = sys.argv[1]
	else:
		e = 1
	q = 0
	if predOff == True:
		b = 0
	if loadFromFile == 'On':
		itime = loadData(data,fileName)
	if plot == 'On':
		if predOff == False:
			dogQuiver, sheepQuiver = initPlot(data, savePlotPng)
		else:
			sheepQuiver = initPlot(data, savePlotPng)
	lastSnapshot = data['t'][itime]
	lastPlot = data['t'][itime]

	while (data['t'][itime] < TF)*(data['alive'][itime].sum() > 2):
		if data['t'][itime]-lastPlot > plotPeriod:
			print data['t'][itime]
			if plot == 'On':
				if predOff == False:
					plotDataPositions(data, itime, dogQuiver, sheepQuiver, savePlotPng)
				else:
					plotDataPositions(data, itime, 'Off', sheepQuiver, savePlotPng)
			lastPlot = data['t'][itime]

		if saveDataH5 == 'On':
			if data['t'][itime]-lastSnapshot > snapshotPeriod:
				print "Saving at "+str(data['t'][itime])
				saveData(data, e, itime)
				lastSnapshot = data['t'][itime]
				print "Wiping old data"
				wipeData(data)
				itime = 4

		if updateMethod == 'Velocity':
			doVelocityStep(data)
		elif updateMethod == 'Acceleration':
			q = doAccelerationStep(data, q)
		else:
			sys.exit('Invalid updateMethod: ' + updateMethod)

		itime+=1

	if saveDataH5 == 'On':
		print "Saving at "+str(data['t'][itime])
		saveData(data, e, itime)
