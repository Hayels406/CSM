import numpy as np
import math
import matplotlib.pyplot as plt
import sys
import h5py

sys.dont_write_bytecode = True
sys.path.insert(0, '.')
from params import *
from scipy.spatial.distance import cdist

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
	data['sheep'] = np.zeros((cachedTimesteps,NP, 2))
	data['sheepVel'] = np.zeros((cachedTimesteps,NP,2))
	data['dogVel'] = np.zeros((cachedTimesteps,2))

	data['dist_id'] = np.zeros((cachedTimesteps,NP))
	data['dist_ij'] = dict()
	data['dist_ij']['min'] = np.zeros(cachedTimesteps)
	data['dist_ij']['max'] = np.zeros(cachedTimesteps)
	data['dist_ij']['mean'] = np.zeros(cachedTimesteps)

	data['t'] = np.zeros(cachedTimesteps)

	if timeStepMethod == 'Adaptive':
		data['q'] = []
	elif timeStepMethod != 'Euler':
		sys.exit("Error: updateMethod not recognized!")

	if wallType == 'Square' :
		data['walls'] = makeSquareWalls(wallTop,wallBottom,wallLeft,wallRight)
		data['dist_iw'] = np.zeros((cachedTimesteps,NP,4))
		data['forceWalls'] = np.zeros((cachedTimesteps,NP,2))
		data['dist_dw'] = np.zeros((cachedTimesteps,4))
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
	data['sheep'][0] = 2*np.array([[x,y] for x in np.arange(0.,np.ceil(np.sqrt(NP)),1.) for y in np.arange(0.,np.ceil(np.sqrt(NP)),1.)])[0:NP]#np.random.rand(NP,2)*3
	#data['sheep'][0][0] = np.array(sheep_init)
	data['sheepVel'][0] = np.array(sheep_vel_init)


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
		sys.exit('Invalid time step method:'+timeStepMethod+' in doVelocitySheep()')

def doAccelerationStep(data):
	#Acceleration of sheep
	distMatrixNotMe = np.array(map(lambda j:data['sheep'][itime,j] - np.delete(data['sheep'][itime], j, 0),range(NP)))
	normStuff = np.transpose(np.tile(np.transpose(np.linalg.norm(distMatrixNotMe, axis = 2)),(2,1,1)))

	preyMinusPred = data['sheep'][itime] - data['dog'][itime]
	normStuff2 = np.transpose(np.tile(np.linalg.norm(preyMinusPred,axis=1),(2,1)))
	data['Gfunc'][itime] = normStuff2**(-1)
	data['sheepAccFlocking'][itime] = 1./NP*((normStuff**(-1) - a*normStuff)*distMatrixNotMe*(normStuff**(-1))).sum(axis=1)
	data['sheepAccFlight'][itime] = data['Gfunc'][itime]*preyMinusPred*(normStuff2**(-1))
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
	#data['dist_ij']['max'][itime] = max(dist_ij)
	#data['dist_ij']['min'][itime] = min(dist_ij)
	#data['dist_ij']['mean'][itime] = np.mean(dist_ij)

	#f_ij = np.zeros((NP, NP, 2))
	#for i in range(NP-1):
	#	f_ij[i, 1+i:, :] = map(lambda j: A*np.exp((sheepSize*2 - dist_ij[i,j])/B)*(data['sheep'][itime][i] - data['sheep'][itime][j])/dist_ij[i,j], range(i+1, NP))

	#data['forceSheep'][itime] = (-np.transpose(f_ij, axes=(1, 0, 2)) + f_ij).sum(axis = 1)
	#print data['forceSheep'][itime]

	#print np.array(map(lambda i:np.array(map(lambda j: A*np.exp((sheepSize*2 - dist_ij[i,j])/B)*(data['sheep'][itime][i] - data['sheep'][itime][j])/dist_ij[i,j], np.delete(range(NP), i, 0))).sum(axis = 0),range(NP)))
	data['forceSheep'][itime] = np.array(map(lambda i: np.nansum(np.fromfunction(lambda j,k: A*np.exp((sheepSize*2 - dist_ij[i,j])/B)*(data['sheep'][itime][i,k] - data['sheep'][itime][j,k])/dist_ij[i,j],(NP,2),dtype=int),axis=0),range(NP)))
	#print data['forceSheep'][itime]
	#print np.shape(np.transpose(data['forceSheep'][itime]))
	#data['forceSheep'][itime] = data['forceSheep'][itime] + np.transpose(data['forceSheep'][itime])

	if len(data['walls'])>0:
		data['sheepAcc'][itime] += data['forceWalls'][itime]/60. + data['forceSheep'][itime]/60.

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



def initPlot(data):
	plt.figure()
	dogTheta = np.zeros(1)
	sheepTheta = np.zeros(1)
	dogQuiver = plt.quiver(data['dog'][0, 0], data['dog'][0, 1], np.cos(dogTheta), np.sin(dogTheta), scale = 30, color = 'red')
	sheepQuiver = plt.quiver(data['sheep'][0,:,0], data['sheep'][0,:,1], np.cos(sheepTheta), np.sin(sheepTheta), scale = 30)
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
		plt.savefig(str(0).zfill(7)+'.png')
	return dogQuiver, sheepQuiver

def plotDataPositions(data):
	dogTheta = np.arctan2(data['dogVel'][itime-1,1], data['dogVel'][itime-1,0])
	sheepTheta = np.arctan2(data['sheepVel'][itime-1,:,1], data['sheepVel'][itime-1,:,0])
	dogQuiver.set_offsets(np.transpose([data['dog'][itime, 0], data['dog'][itime, 1]]))
	dogQuiver.set_UVC(np.cos(dogTheta),np.sin(dogTheta))
	sheepQuiver.set_offsets(np.transpose([data['sheep'][itime,:, 0], data['sheep'][itime,:, 1]]))
	sheepQuiver.set_UVC(np.cos(sheepTheta), np.sin(sheepTheta))
	plt.pause(0.005)
	if savePlotPng == 'On':
		plt.savefig(str(itime).zfill(7)+'.png')


def saveData(data):
	h5f = h5py.File('data-%07d.h5'%itime, 'w')
	for k in data.keys():
		if type(data[k]) is dict:
			for k2 in data[k].keys():
				try:
					h5f.create_dataset(k+'-'+k2, data=data[k][k2],compression="gzip")
				except:
					h5f.create_dataset(k+'-'+k2, data=data[k][k2])
		else:
			try:
				h5f.create_dataset(k, data=data[k],compression="gzip")
			except:
				h5f.create_dataset(k, data=data[k])
	h5f.close()

def loadData(data,fn):
	h5f = h5py.File(fn,'r')
	itime = int(fn[5:-3])
	data['dog'] = np.copy(h5f['dog'])
	data['sheep'] = np.copy(h5f['sheep'])
	data['sheepVel'] = np.copy(h5f['sheepVel'])
	data['dogVel'] = np.copy(h5f['dogVel'])

	data['dist_id'] = np.copy(h5f['dist_id'])
	data['dist_ij'] = dict()
	data['dist_ij']['min'] = np.copy(h5f['dist_ij-min'])
	data['dist_ij']['max'] = np.copy(h5f['dist_ij-max'])
	data['dist_ij']['mean'] = np.copy(h5f['dist_ij-mean'])

	data['t'] = np.copy(h5f['t'])

	if timeStepMethod == 'Adaptive':
		data['q'] = np.copy(h5f['q']).tolist()
	elif timeStepMethod != 'Euler':
		sys.exit("Error: updateMethod not recognized!")

	if wallType == 'Square' :
		data['walls']['eqn'] = np.copy(h5f['walls-eqn'])
		data['walls']['n'] = np.copy(h5f['walls-n'])
		data['dist_iw'] = np.copy(h5f['dist_iw'])
		data['forceWalls'] = np.copy(h5f['forceWalls'])
		data['dist_dw'] = np.copy(h5f['dist_dw'])
	else:
		data['walls'] = []

	if updateMethod == 'Velocity':
		data['bterm'] = np.copy(h5f['bterm'])
		data['aterm'] = np.copy(h5f['aterm'])
	elif updateMethod == 'Acceleration':
		data['sheepAccFlocking'] = np.copy(h5f['sheepAccFlocking'])
		data['sheepAccFlight'] = np.copy(h5f['sheepAccFlight'])
		data['sheepAcc'] = np.copy(h5f['sheepAcc'])
		data['dogAcc'] = np.copy(h5f['dogAcc'])
		data['Gfunc'] = np.copy(h5f['Gfunc'])
		data['Hfunc'] = np.copy(h5f['Hfunc'])
		data['Ffunc'] = np.copy(h5f['Ffunc'])
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
					data[k][k2][5:] = data[k][k2][5:]*0
			else:
				data[k][0:5] = data[k][itime-4:itime+1]
				data[k][5:] = data[k][5:]*0


###############
#Simulation
###############
plt.close()
data = init()
itime = 0
initCond(data)
if loadFromFile == 'On':
	itime = loadData(data,fileName)
dogQuiver, sheepQuiver = initPlot(data)
lastSnapshot = data['t'][itime]
lastPlot = data['t'][itime]

while data['t'][itime] < TF:
	if data['t'][itime]-lastPlot > plotPeriod:
		print data['t'][itime]
		print itime
		plotDataPositions(data)
		lastPlot = data['t'][itime]

	if data['t'][itime]-lastSnapshot > snapshotPeriod:
		saveData(data)
		print "Saving at "+str(data['t'][itime])
		lastSnapshot = data['t'][itime]
		wipeData(data)
		itime = 4

	if updateMethod == 'Velocity':
		doVelocityStep(data)
	elif updateMethod == 'Acceleration':
		doAccelerationStep(data)
	else:
		sys.exit('Invalid updateMethod: '+updateMethod)

	itime+=1

#save as matlab .mat file
#import scipy.io
#
#x = np.linspace(0, 2 * np.pi, 100)
#y = np.cos(x)
#
#scipy.io.savemat('test.mat', dict(x=x, y=y))
