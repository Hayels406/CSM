import numpy as np
import math
import matplotlib.pyplot as plt
import sys
import h5py
import os
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial import Voronoi
from scipy.spatial.distance import cdist
sys.dont_write_bytecode = True
sys.path.insert(0,os.getcwd())
from params import *



def makeSquareWalls(wallTop,wallBottom,wallLeft,wallRight):
	#Walls defined by four lines: Ax+By+C = 0
	w = dict()
	w['eqn'] = [[0,1,-wallTop],[0,1,-wallBottom],[1,0,-wallLeft],[1,0,-wallRight]]
	w['n'] = [np.array([0,-1]),np.array([0,1]),np.array([1,0]),np.array([-1,0])]
	return w

def init():
	#Sets up main data structure
	#Returns: dictionary of simulation data
	data = dict()
	cachedTimesteps = 5*int(round(TF/dt))
	data['dog'] = np.zeros((cachedTimesteps, 2))
	data['dogVel'] = np.zeros((cachedTimesteps,2))
	data['sheep'] = np.zeros((cachedTimesteps,NP, 2))
	data['sheepVel'] = np.zeros((cachedTimesteps,NP,2))

	data['dist_id'] = np.zeros((cachedTimesteps,NP))
	data['dist_ij'] = dict()
	data['dist_ij']['min'] = np.zeros(cachedTimesteps)
	data['dist_ij']['max'] = np.zeros(cachedTimesteps)
	data['dist_ij']['mean'] = np.zeros(cachedTimesteps)

	data['t'] = np.zeros(cachedTimesteps)
	data['dog_index_neighbours'] = []

	if timeStepMethod == 'Adaptive':
		data['q'] = []
	elif timeStepMethod != 'Euler':
		sys.exit("Error: updateMethod not recognized!")

	if wallType == 'Square' :
		data['walls'] = makeSquareWalls(wallTop,wallBottom,wallLeft,wallRight)
		data['dist_iw'] = np.zeros((cachedTimesteps,NP,4))
		data['forceWalls'] = np.zeros((cachedTimesteps,NP,2))
		data['dist_dw'] = np.zeros((cachedTimesteps,4))
		if sheepSheepInteration == 'On':
			data['forceSheep'] = np.zeros((cachedTimesteps, NP, 2))
	else:
		data['walls'] = []

	if updateMethod == 'Velocity':
		data['bterm'] = np.zeros((cachedTimesteps,NP,2))
		data['aterm'] = np.zeros((cachedTimesteps,NP,2))
	elif updateMethod == 'Acceleration':
		data['sheepAccFlocking'] = np.zeros((cachedTimesteps,NP,2))
		data['sheepAccFlight'] = np.zeros((cachedTimesteps,NP,2))
		data['sheepAcc'] = np.zeros((cachedTimesteps,NP,2))
		data['dogAcc'] = np.zeros((cachedTimesteps,2))
		data['Gfunc'] = np.zeros((cachedTimesteps,NP,2))
		data['Hfunc'] = np.zeros((cachedTimesteps,NP,2))
		data['Ffunc'] = np.zeros((cachedTimesteps,NP,2))
	else:
		sys.exit("Error: updateMethod not recognized!")

	return data

def initCond(data):
	data['dog'][0] = np.array(dog_init)
	data['dogVel'][0] = np.array(dog_vel_init)
	if sheep_init == 'Grid':
		data['sheep'][0] = 2*np.array([[x,y] for x in np.arange(0.,np.ceil(np.sqrt(NP)),1.) for y in np.arange(0.,np.ceil(np.sqrt(NP)),1.)])[0:NP]#np.random.rand(NP,2)*3
	elif sheep_init == 'Random':
		min_dist_ij = 0.
		while min_dist_ij < 1.0:
			data['sheep'][0,:,0] = np.random.rand(NP)*(wallRight - wallLeft - 5.) + wallLeft
			data['sheep'][0,:,1] = np.random.rand(NP)*(wallTop - wallBottom/2. - 5.) + wallBottom/2.
			dist_ij = cdist(data['sheep'][itime], data['sheep'][itime])
			min_dist_ij = np.min(dist_ij[np.nonzero(dist_ij)])
			print min_dist_ij
	sheepTheta = np.random.rand(NP)*2*math.pi
	data['sheepVel'][0,:,0] = sheep_vel_init*np.cos(sheepTheta)
	data['sheepVel'][0,:,1] = sheep_vel_init*np.sin(sheepTheta)


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

	data['dist_id'][itime] = cdist(data['sheep'][itime], data['dog'][itime].reshape(1,2)).reshape(NP)
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
			data['q'].append(np.floor(100.*q_f)/100.)
			data['t'][itime+1] = data['t'][itime] + dt
		else:
			sheepTemp = data['sheep'][itime] + data['q'][-1]*(dt/24.)*(55.*data['sheepVel'][itime] - 59.*data['sheepVel'][itime - 1] + 37.*data['sheepVel'][itime - 2] - 9.*data['sheepVel'][itime - 3])
			sheepVelTemp = (sheepTemp - data['sheep'][itime])/(dt*data['q'][-1])
			data['sheep'][itime + 1] = data['sheep'][itime] + data['q'][-1]*(dt/24.)*(9.*sheepVelTemp + 19.*data['sheepVel'][itime] - 5.*data['sheepVel'][itime - 1] + data['sheepVel'][itime - 2])
			data['dog'][itime + 1] = data['dog'][itime] + data['dogVel'][itime]*dt*data['q'][-1]
			q_f = np.min(((270./19.)*((epsilon*dt)/(np.abs(data['sheep'][itime + 1] - sheepTemp))))**(0.25))
			data['t'][itime+1] = data['t'][itime] + dt*data['q'][-1]
			data['q'].append(np.floor(100.*q_f)/100.)

	else:
		sys.exit('Invalid time step method:'+timeStepMethod+' in doVelocityStep()')

def doAccelerationStep(data):
	#Acceleration of sheep
	distMatrixNotMe = np.array(map(lambda j:data['sheep'][itime,j] - np.delete(data['sheep'][itime], j, 0),range(NP)))
	normStuff = np.transpose(np.tile(np.transpose(np.linalg.norm(distMatrixNotMe, axis = 2)),(2,1,1)))

	preyMinusPred = data['sheep'][itime] - data['dog'][itime]
	normStuff2 = np.transpose(np.tile(np.linalg.norm(preyMinusPred,axis=1),(2,1)))
	data['Gfunc'][itime] = normStuff2**(-1)
	data['sheepAccFlocking'][itime] = f/NP*((normStuff**(-1) - a*normStuff)*distMatrixNotMe*(normStuff**(-1))).sum(axis=1)

	if allSheepSeeDog == 'On':
		data['sheepAccFlight'][itime] = data['Gfunc'][itime]*preyMinusPred*(normStuff2**(-1))
	elif allSheepSeeDog == 'Off':
		vor = Voronoi(np.concatenate((data['sheep'][itime], data['dog'][itime][None,:]), axis = 0))
		data['dog_index_neighbours'] = np.array([t[1] for t in [(b1, a1) for a1, b1 in vor.ridge_dict.keys()] + vor.ridge_dict.keys() if t[0] == NP]) #Calculates the Voronoi neighbours of dog
		data['sheepAccFlight'][itime][data['dog_index_neighbours']] = (data['Gfunc'][itime]*preyMinusPred*(normStuff2**(-1)))[data['dog_index_neighbours']]
	else:
		sys.exit('Invalid allSheepSeeDog:'+allSheepSeeDog+' in doAccelerationStep()')

	data['sheepAcc'][itime] = -data['sheepVel'][itime] + data['sheepAccFlight'][itime]+ data['sheepAccFlocking'][itime]

	#Acceleration of dog
	predMinusPrey = -(data['dog'][itime] - data['sheep'][itime])
	normStuff3 = np.transpose(np.tile(np.linalg.norm(predMinusPrey,axis=1),(2,1)))
	data['Hfunc'][itime] = 1./(normStuff3**1. + 1.)
	data['dogAcc'][itime] = -data['dogVel'][itime] + 1./NP*(predMinusPrey*(normStuff3**(-1))*data['Hfunc'][itime]).sum(axis=0)

	for wall in range(len(data['walls']['eqn'])):
		w = data['walls']['eqn'][wall]
		data['dist_iw'][itime,:,wall] = np.abs(w[0]*data['sheep'][itime,:,0] + w[1]*data['sheep'][itime,:,1] + w[2])/np.sqrt(w[0]**2 + w[1]**2)#np.sqrt(((sheep[itime] - wall_point)**2).sum(axis = 0))
		data['forceWalls'][itime,:,0] += A*np.exp((sheepSize - data['dist_iw'][itime,:,wall])/B)*data['walls']['n'][wall][0]
		data['forceWalls'][itime,:,1] += A*np.exp((sheepSize - data['dist_iw'][itime,:,wall])/B)*data['walls']['n'][wall][1]

	dist_ij = cdist(data['sheep'][itime], data['sheep'][itime])
	data['dist_ij']['max'][itime] = np.max(dist_ij)
	data['dist_ij']['min'][itime] = np.min(dist_ij[np.nonzero(dist_ij)])
	data['dist_ij']['mean'][itime] = np.mean(dist_ij)
	print data['dist_ij']['min'][itime]

	data['dist_id'][itime] = cdist(data['sheep'][itime], data['dog'][itime].reshape(1,2)).reshape(NP)
	if sheepSheepInteration == 'On':
		data['forceSheep'][itime] = np.array(map(lambda i: np.nansum(np.fromfunction(lambda j,k: A*np.exp((sheepSize*2 - dist_ij[i,j])/B)*(data['sheep'][itime][i,k] - data['sheep'][itime][j,k])/dist_ij[i,j],(NP,2),dtype=int),axis=0),range(NP)))

	if len(data['walls'])>0:
		if sheepSheepInteration == 'On':
			data['sheepAcc'][itime] += data['forceWalls'][itime]/60. + data['forceSheep'][itime]/60.
		else:
			data['sheepAcc'][itime] += data['forceWalls'][itime]/60.

	if timeStepMethod == 'Euler':
		#Euler time step
		data['dogVel'][itime + 1] = data['dogVel'][itime] + data['dogAcc'][itime]*dt
		data['sheepVel'][itime + 1] = data['sheepVel'][itime] + data['sheepAcc'][itime]*dt
		data['dog'][itime + 1] = data['dog'][itime] + data['dogVel'][itime]*dt
		data['sheep'][itime + 1] = data['sheep'][itime] + data['sheepVel'][itime]*dt
		data['t'][itime+1] = data['t'][itime] + dt
	elif timeStepMethod == 'Adaptive':
		#Adaptive time step
		if itime == 0:
			data['dogVel'][itime + 1] = data['dogVel'][itime] + data['dogAcc'][itime]*dt
			data['sheepVel'][itime + 1] = data['sheepVel'][itime] + data['sheepAcc'][itime]*dt
			data['sheep'][itime + 1] = data['sheep'][itime] + data['sheepVel'][itime]*dt
			data['dog'][itime + 1] = data['dog'][itime] + data['dogVel'][itime]*dt
			data['t'][itime+1] = data['t'][itime] + dt
		elif itime == 1:
			data['sheepVel'][itime + 1] = data['sheepVel'][itime] + 1.5*dt*data['sheepAcc'][itime] - 0.5*dt*data['sheepAcc'][itime - 1]
			data['dogVel'][itime +1] = data['dogVel'][itime] + data['dogAcc'][itime]*dt
			data['sheep'][itime + 1] = data['sheep'][itime] + data['sheepVel'][itime]*dt
			data['dog'][itime + 1] = data['dog'][itime] + data['dogVel'][itime]*dt
			data['t'][itime+1] = data['t'][itime] + dt
		elif itime == 2:
			data['sheepVel'][itime + 1] = data['sheepVel'][itime] + dt*((23./12.)*data['sheepAcc'][itime] - (4./3.)*data['sheepAcc'][itime - 1] + (5./12.)*data['sheepAcc'][itime - 2])
			data['dogVel'][itime + 1] = data['dogVel'][itime] + data['dogAcc'][itime]*dt
			data['sheep'][itime + 1] = data['sheep'][itime] + data['sheepVel'][itime]*dt
			data['dog'][itime + 1] = data['dog'][itime] + data['dogVel'][itime]*dt
			data['t'][itime+1] = data['t'][itime] + dt
		elif itime == 3:
			sheepVelTemp = data['sheepVel'][itime] + (dt/24.)*(55.*data['sheepAcc'][itime] - 59.*data['sheepAcc'][itime - 1] + 37.*data['sheepAcc'][itime - 2] - 9.*data['sheepAcc'][itime - 3])
			sheepAccTemp = (sheepVelTemp - data['sheepVel'][itime])/dt
			data['sheepVel'][itime + 1] = data['sheepVel'][itime] + (dt/24.)*(9.*sheepAccTemp + 19.*data['sheepAcc'][itime] - 5.*data['sheepAcc'][itime - 1] + data['sheepAcc'][itime - 2])
			data['dogVel'][itime + 1] = data['dogVel'][itime] + data['dogAcc'][itime]*dt
			data['sheep'][itime + 1] = data['sheep'][itime] + data['sheepVel'][itime]*dt
			data['dog'][itime + 1] = data['dog'][itime] + data['dogVel'][itime]*dt
			q_f = np.min(((270./19.)*((epsilon*dt)/(np.abs(data['sheepVel'][itime + 1] - sheepVelTemp))))**(0.25))
			data['q'].append(np.floor(100.*q_f)/100.)
			data['t'][itime+1] = data['t'][itime] + dt
		else:
			sheepVelTemp = data['sheepVel'][itime] + data['q'][-1]*(dt/24.)*(55.*data['sheepAcc'][itime] - 59.*data['sheepAcc'][itime - 1] + 37.*data['sheepAcc'][itime - 2] - 9.*data['sheepAcc'][itime - 3])
			sheepAccTemp = (sheepVelTemp - data['sheepVel'][itime])/(dt*data['q'][-1])
			data['sheepVel'][itime + 1] = data['sheepVel'][itime] + data['q'][-1]*(dt/24.)*(9.*sheepAccTemp + 19.*data['sheepAcc'][itime] - 5.*data['sheepAcc'][itime - 1] + data['sheepAcc'][itime - 2])
			data['dogVel'][itime + 1] = data['dogVel'][itime] + data['dogAcc'][itime]*dt*data['q'][-1]
			data['sheep'][itime + 1] = data['sheep'][itime] + data['sheepVel'][itime]*dt*data['q'][-1]
			data['dog'][itime + 1] = data['dog'][itime] + data['dogVel'][itime]*dt*data['q'][-1]
			q_f = np.min(((270./19.)*((epsilon*dt)/(np.abs(data['sheepVel'][itime + 1] - sheepVelTemp))))**(0.25))
			data['t'][itime+1] = data['t'][itime] + dt*data['q'][-1]
			data['q'].append(np.floor(100.*q_f)/100.)
	else:
		sys.exit('Invalid time step method:'+timeStepMethod+' in doAccelerationStep()')



def initPlot(data, savePlotPng):
	plt.figure()
	dogTheta = np.zeros(1)
	sheepTheta = np.zeros(1)
	from matplotlib import colors
	cmap = colors.ListedColormap(['black', 'blue'])
	bounds=[0,0.5,1.0]
	norm = colors.BoundaryNorm(bounds, cmap.N)
	dogQuiver = plt.quiver(data['dog'][0, 0], data['dog'][0, 1], np.cos(dogTheta), np.sin(dogTheta), scale = 30, color = 'red')
	sheepQuiver = plt.quiver(data['sheep'][0,:,0], data['sheep'][0,:,1], np.cos(sheepTheta), np.sin(sheepTheta), scale = 30, cmap = cmap, norm = norm)
	plt.axis([wallLeft,wallRight,wallBottom,wallTop])
	plt.axes().set_aspect('equal')
	plotid = 0
	if wallType == 'Square':
		plt.axhline(wallTop, color = 'r', lw = 5)
		plt.axhline(wallBottom, color = 'r', lw = 5)
		plt.axvline(wallLeft, color = 'r', lw = 5)
		plt.axvline(wallRight, color = 'r', lw = 5)
	plt.pause(0.005)
	if savePlotPng == 'On':
		plt.savefig('frames/'+str(0).zfill(7)+'.png')
	return dogQuiver, sheepQuiver

def plotDataPositions(data, tstep, dQuiv, sQuiv, savePlotPng):
	colorQuiver = np.zeros(NP)
	if showSheepDogCanSee == 'On':
		colorQuiver[data['dog_index_neighbours']] = 1
	dogTheta = np.arctan2(data['dogVel'][tstep-1,1], data['dogVel'][tstep-1,0])
	sheepTheta = np.arctan2(data['sheepVel'][tstep-1,:,1], data['sheepVel'][tstep-1,:,0])
	dQuiv.set_offsets(np.transpose([data['dog'][tstep, 0], data['dog'][tstep, 1]]))
	dQuiv.set_UVC(np.cos(dogTheta),np.sin(dogTheta))
	sQuiv.set_offsets(np.transpose([data['sheep'][tstep,:, 0], data['sheep'][tstep,:, 1]]))
	sQuiv.set_UVC(np.cos(sheepTheta), np.sin(sheepTheta), colorQuiver)
	plt.pause(0.005)
	if savePlotPng == 'On':
		plt.savefig('frames/'+str(np.floor(data['t'][tstep])).zfill(7)+'.png')


def saveData(data):
	h5f = h5py.File('data-%07d.h5'%int(data['t'][itime]), 'w')
	h5f.create_dataset('itime', data = [itime])
	for k in data.keys():
		if type(data[k]) is dict:
			for k2 in data[k].keys():
				try:
					h5f.create_dataset(k+'-'+k2, data=data[k][k2][0:itime+1],compression="gzip")
				except:
					h5f.create_dataset(k+'-'+k2, data=data[k][k2][0:itime+1])
		else:
			try:
				h5f.create_dataset(k, data=data[k][0:itime+1],compression="gzip")
			except:
				h5f.create_dataset(k, data=data[k][0:itime+1])
	h5f.close()

def loadData(data,fn):
	h5f = h5py.File(fn,'r')
	itime = np.copy(h5f['itime'])[0]
	data['dog'][0:itime+1] = np.copy(h5f['dog'])
	data['sheep'][0:itime+1] = np.copy(h5f['sheep'])
	data['sheepVel'][0:itime+1] = np.copy(h5f['sheepVel'])
	data['dogVel'][0:itime+1] = np.copy(h5f['dogVel'])

	data['dist_id'][0:itime+1] = np.copy(h5f['dist_id'])
	data['dist_ij']['min'][0:itime+1] = np.copy(h5f['dist_ij-min'])
	data['dist_ij']['max'][0:itime+1] = np.copy(h5f['dist_ij-max'])
	data['dist_ij']['mean'][0:itime+1] = np.copy(h5f['dist_ij-mean'])

	data['t'][0:itime+1] = np.copy(h5f['t'])

	if timeStepMethod == 'Adaptive':
		data['q'] = np.copy(h5f['q']).tolist()
	elif timeStepMethod != 'Euler':
		sys.exit("Error: updateMethod not recognized!")

	if wallType == 'Square' :
		data['walls']['eqn'] = np.copy(h5f['walls-eqn'])
		data['walls']['n'] = np.copy(h5f['walls-n'])
		data['dist_iw'][0:itime+1] = np.copy(h5f['dist_iw'])
		data['forceWalls'][0:itime+1] = np.copy(h5f['forceWalls'])
		data['dist_dw'][0:itime+1] = np.copy(h5f['dist_dw'])
	else:
		data['walls'] = []

	if updateMethod == 'Velocity':
		data['bterm'][0:itime+1] = np.copy(h5f['bterm'])
		data['aterm'][0:itime+1] = np.copy(h5f['aterm'])
	elif updateMethod == 'Acceleration':
		data['sheepAccFlocking'][0:itime+1] = np.copy(h5f['sheepAccFlocking'])
		data['sheepAccFlight'][0:itime+1] = np.copy(h5f['sheepAccFlight'])
		data['sheepAcc'][0:itime+1] = np.copy(h5f['sheepAcc'])
		data['dogAcc'][0:itime+1] = np.copy(h5f['dogAcc'])
		data['Gfunc'][0:itime+1] = np.copy(h5f['Gfunc'])
		data['Hfunc'][0:itime+1] = np.copy(h5f['Hfunc'])
		data['Ffunc'][0:itime+1] = np.copy(h5f['Ffunc'])
	else:
		sys.exit("Error: updateMethod not recognized!")
	return itime

def wipeData(data):
	for k in data.keys():
		if k != 'walls':
			if k=='q':
				data['q'] = [data['q'][-1]]
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
	if loadFromFile == 'On':
		itime = loadData(data,fileName)
	if plot == 'On':
		dogQuiver, sheepQuiver = initPlot(data, savePlotPng)
	lastSnapshot = data['t'][itime]
	lastPlot = data['t'][itime]

	while data['t'][itime] < TF:
		if data['t'][itime]-lastPlot > plotPeriod:
			print data['t'][itime]
			if plot == 'On':
				plotDataPositions(data, itime, dogQuiver, sheepQuiver, savePlotPng)
			lastPlot = data['t'][itime]

		if saveDataH5 == 'On':
			if data['t'][itime]-lastSnapshot > snapshotPeriod:
				print "Saving at "+str(data['t'][itime])
				saveData(data)
				lastSnapshot = data['t'][itime]
				print "Wiping old data"
				wipeData(data)
				itime = 4

		if updateMethod == 'Velocity':
			doVelocityStep(data)
		elif updateMethod == 'Acceleration':
			doAccelerationStep(data)
		else:
			sys.exit('Invalid updateMethod: ' + updateMethod)

		itime+=1

	if saveDataH5 == 'On':
		print "Saving at "+str(data['t'][itime])
		saveData(data)
