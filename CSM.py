import numpy as np
import math
import matplotlib.pyplot as plt
import sys
import json
sys.dont_write_bytecode = True
sys.path.insert(0, '')
from params import *

def makeSquareWalls(wallTop,wallBottom,wallLeft,wallRight):
	#Walls defined by four lines: Ax+By+C = 0
	#Returns: List of 4 lists
	return [[0,1,-wallTop],[0,1,-wallBottom],[1,0,-wallLeft],[1,0,-wallRight]]

def init():
	#Sets up main data structure
	#Returns: dictionary of simulation data
	data = dict()
	data['cachedTimesteps'] = 10*int(round(TF/dt))
	data['dog'] = np.zeros((data['cachedTimesteps'], 2))
	data['sheep'] = np.zeros((data['cachedTimesteps'], 2))
	data['sheepVel'] = np.zeros((data['cachedTimesteps'],2))
	data['dogVel'] = np.zeros((data['cachedTimesteps'],2))

	data['dist_id'] = np.zeros((data['cachedTimesteps'],NP))
	data['dist_ij'] = dict()
	data['dist_ij']['min'] = np.zeros(data['cachedTimesteps'])
	data['dist_ij']['max'] = np.zeros(data['cachedTimesteps'])
	data['dist_ij']['mean'] = np.zeros(data['cachedTimesteps'])

	data['t'] = np.zeros(data['cachedTimesteps'])
	
	if timeStepMethod == 'Adaptive':
		data['q'] = []
	elif timeStepMethod != 'Euler':
		sys.exit("Error: updateMethod not recognized!")

	if wallType == 'Square' :
		data['walls'] = makeSquareWalls(wallTop,wallBottom,wallLeft,wallRight)
		data['dist_iw'] = np.zeros((data['cachedTimesteps'],NP,4))
		data['forceWalls'] = np.zeros((data['cachedTimesteps'],2))
		data['dist_dw'] = np.zeros((data['cachedTimesteps'],4))
	else:
		data['walls'] = []

	if updateMethod == 'Velocity':
		data['bterm'] = np.zeros((data['cachedTimesteps'],2))
		data['aterm'] = np.zeros(data['cachedTimesteps'])
	elif updateMethod == 'Acceleration':
		data['sheepAccFlocking'] = np.zeros((data['cachedTimesteps'],2))
		data['sheepAccFlight'] = np.zeros((data['cachedTimesteps'],2))
		data['sheepAcc'] = np.zeros((data['cachedTimesteps'],2))
		data['dogAcc'] = np.zeros((data['cachedTimesteps'],2))
		data['Gfunc'] = np.zeros((data['cachedTimesteps'],2))
		data['Hfunc'] = np.zeros((data['cachedTimesteps'],2))
	else:
		sys.exit("Error: updateMethod not recognized!")

	
	return data

def initCond(data):
	data['dog'][0] = np.array(dog_init)
	data['sheep'][0] = np.array(sheep_init)
	data['sheepVel'][0] = np.array(sheep_vel_init)
	data['dogVel'][0] = np.array(dog_vel_init)


###############
#Velocity
###############
def doVelocityStep(data):
	#Velocity of sheep
    #distMatrixNotMe = np.array(map(lambda j:sheep[itime,j] - np.delete(sheep[itime], j, 0),range(NP)))
    #normStuff[itime] = np.transpose(np.tile(np.transpose(np.linalg.norm(distMatrixNotMe, axis = 2)**(-2)),(2,1,1)))
    preyMinusPred = data['sheep'][itime] - data['dog'][itime]
    normstuff2 = np.transpose(np.tile(np.linalg.norm(preyMinusPred),(2,1)))
    data['bterm'][itime] = b*preyMinusPred*normstuff2**(-2)
    data['sheepVel'][itime] = data['bterm'][itime]#(distMatrixNotMe*(normStuff2 - a)).sum(axis = 1)/NP + bterm[itime]

    #Velocity of dog
    distDogSheep = data['sheep'][itime] - data['dog'][itime]
    frac = distDogSheep*np.transpose(np.tile((np.linalg.norm(distDogSheep))**(-p),(2,1)))
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
    preyMinusPred = data['sheep'][itime] - data['dog'][itime]
    normStuff2 = np.transpose(np.tile(np.linalg.norm(preyMinusPred),(2,1)))
    data['Gfunc'][itime] = normStuff2**(-1)
    data['sheepAccFlight'][itime] = data['Gfunc'][itime]*preyMinusPred*(normStuff2**(-1))
    data['sheepAcc'][itime] = -data['sheepVel'][itime] + data['sheepAccFlight'][itime]
    
    #Acceleration of dog
    predMinusPrey = -(data['dog'][itime] - data['sheep'][itime])
    normStuff3 = np.transpose(np.tile(np.linalg.norm(predMinusPrey),(2,1)))
    data['Hfunc'][itime] = 1./(normStuff3**1. + 1.)
    data['dogAcc'][itime] = -data['dogVel'][itime] + predMinusPrey*(normStuff3**(-1))*data['Hfunc'][itime]

    for wall in range(len(data['walls'])):
    	w = data['walls'][wall]
        wall_point = [0.,0.]
        wall_point[0] = (w[0]*w[2] - data['sheep'][itime,1]*w[0]*w[1] + data['sheep'][itime,0]*w[1]**2)/(w[0]**2 + w[1]**2)
        wall_point[1] = -(w[0]/(w[1]+0.01))*wall_point[0] - w[2]/(w[1]+0.01)
        #wall_point = [(y*(y*sheep[itime, 0] - x*sheep[itime, 1]) + x*z)/(x**2 + y**2), (-z*(y*sheep[itime, 0] - x*sheep[itime, 1]) - y*z)/(x**2 + y**2)]
        data['dist_iw'][itime][0][wall] = abs(w[0]*data['sheep'][itime,0] + w[1]*data['sheep'][itime,1] + w[2])/np.sqrt(w[0]**2 + w[1]**2)#np.sqrt(((sheep[itime] - wall_point)**2).sum(axis = 0))
        n_iw = (data['sheep'][itime] - wall_point)/data['dist_iw'][itime][0][wall]
        data['forceWalls'][itime] += A*np.exp((sheepSize - data['dist_iw'][itime][0][wall])/B)*n_iw
    
    if len(data['walls'])>0:
    	data['sheepAcc'][itime] += data['forceWalls'][itime] 

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
	sheepQuiver = plt.quiver(data['sheep'][0, 0], data['sheep'][0, 1], np.cos(sheepTheta), np.sin(sheepTheta), scale = 30)
	plt.axis([wallLeft,wallRight,wallBottom,wallTop])
	plt.axes().set_aspect('equal')
	plotid = 0
	if wallType == 'Square':
		plt.axhline(wallTop, color = 'r', lw = 5)
		plt.axhline(wallBottom, color = 'r', lw = 5)
		plt.axvline(wallLeft, color = 'r', lw = 5)
		plt.axvline(wallRight, color = 'r', lw = 5)
	plt.pause(0.005)
	return dogQuiver, sheepQuiver

def plotDataPositions(data):
	dogTheta = np.arctan2(data['dogVel'][itime,1], data['dogVel'][itime,0])
	sheepTheta = np.arctan2(data['sheepVel'][itime,1], data['sheepVel'][itime,0])
	dogQuiver.set_offsets(np.transpose([data['dog'][itime, 0], data['dog'][itime, 1]]))
	dogQuiver.set_UVC(np.cos(dogTheta),np.sin(dogTheta))
	sheepQuiver.set_offsets(np.transpose([data['sheep'][itime, 0], data['sheep'][itime, 1]]))
	sheepQuiver.set_UVC(np.cos(sheepTheta), np.sin(sheepTheta))
	plt.pause(0.005)

###############
#Simulation
###############
plt.close()
data = init()
initCond(data)
dogQuiver, sheepQuiver = initPlot(data)
itime = 0

while data['t'][itime] < TF:
	if itime % 1000 == 0:
		print data['t'][itime]
		print data['sheep'][itime]
		plotDataPositions(data)

	if itime % 10000 == 0:
		np.save('%d.npy'%itime,data)

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













