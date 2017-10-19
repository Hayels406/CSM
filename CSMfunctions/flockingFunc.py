import numpy as np
import sys
from scipy.stats import itemfreq
from sklearn.neighbors import KDTree

def doPredatorPrey(DATA, ITIME, Np, F, A, GROUPSIZE, SHEEPSIZE, PREDATION, GAUSSIAN):
    tree = KDTree(DATA['sheep'][ITIME][DATA['alive'][ITIME]])
    idx = tree.query(DATA['sheep'][ITIME][DATA['alive'][ITIME]], np.min([GROUPSIZE, sum(DATA['alive'][ITIME])]))[1][:,1:]
    distMatrix = np.array(map(lambda j:DATA['sheep'][ITIME,j] - DATA['sheep'][ITIME][DATA['alive'][ITIME]],np.array(range(Np))[DATA['alive'][ITIME]]))
    if PREDATION == 'Off':
        normStuff = np.linalg.norm(distMatrix, axis = 2).reshape(Np,Np,1) + SHEEPSIZE
    else:
        normStuff = np.linalg.norm(distMatrix, axis = 2).reshape(sum(DATA['alive'][ITIME]), sum(DATA['alive'][ITIME]), 1) + SHEEPSIZE
    if GAUSSIAN == 'On':
        sheepAccFlocking = (1./np.min([GROUPSIZE, sum(DATA['alive'][ITIME])]))*np.array(map(lambda j:((F*normStuff[j,idx[j],:]**(-1) - A*normStuff[j,idx[j],:])*np.exp(-normStuff[j,idx[j],:]**2/VISUALDIST**2)*distMatrix[j,idx[j],:]*(normStuff[j,idx[j],:]**(-1))), range(sum(DATA['alive'][ITIME])))).sum(axis = 1)
    else:
        sheepAccFlocking = (1./np.min([GROUPSIZE, sum(DATA['alive'][ITIME])]))*np.array(map(lambda j:((F*normStuff[j,idx[j],:]**(-1) - A*normStuff[j,idx[j],:])*distMatrix[j,idx[j],:]*(normStuff[j,idx[j],:]**(-1))), range(sum(DATA['alive'][ITIME])))).sum(axis = 1)
    return sheepAccFlocking




def doPPPoisson(DATA, ITIME, Np, F, A, LAM, GROUPSIZE, PREDATION):
    distributionN = GROUPSIZE + np.random.poisson(LAM, Np) - np.random.poisson(LAM, Np)
    while any(distributionN < 0):
        distributionN = GROUPSIZE + np.random.poisson(LAM, Np) - np.random.poisson(LAM, Np)
    neighCalcTree = itemfreq(distributionN)[:,0]
    neighCalcTree = [int(x) for x in neighCalcTree]
    tree = KDTree(DATA['sheep'][ITIME][DATA['alive'][ITIME]])
    idx = tree.query(DATA['sheep'][ITIME][DATA['alive'][ITIME]], np.min([np.max(neighCalcTree), sum(DATA['alive'][ITIME])]))[1][:,1:]
    distMatrix = np.array(map(lambda j:DATA['sheep'][ITIME,j] - DATA['sheep'][ITIME][DATA['alive'][ITIME]],np.array(range(Np))[DATA['alive'][ITIME]]))
    if PREDATION == 'Off':
        normStuff = np.linalg.norm(distMatrix, axis = 2).reshape(Np,Np,1) + SHEEPSIZE
    else:
        normStuff = np.linalg.norm(distMatrix, axis = 2).reshape(sum(DATA['alive'][ITIME]), sum(DATA['alive'][ITIME]), 1) + SHEEPSIZE

    acc_temp = np.zeros(Np)
    for neighbours in neighCalcTree:
        if neighbours == 0:
            acc_temp[distributionN == neighbours] = 0.0
        else:
			acc_temp[distributionN == neighbours] = (1./np.min([neighbours, sum(DATA['alive'][ITIME])]))*np.array(map(lambda j:((F*normStuff[j,idx[j][:neighbours],:]**(-1) - A*normStuff[j,idx[j][:neighbours],:] )*distMatrix[j,idx[j][:neighbours],:]*(normStuff[j,idx[j][:neighbours],:]**(-1))), range(sum(DATA['alive'][ITIME])))).sum(axis = 1)
    sheepAccFlocking = acc_temp[DATA['alive'][ITIME]]
    return sheepAccFlocking

def doVicsek(DATA, ITIME, Np, R):
    sheepAccFlocking = np.zeros(([DATA['alive'][ITIME]].sum(), 2))
    tree = KDTree(DATA['sheep'][ITIME][DATA['alive'][ITIME]])
    idx = tree.query_radius(DATA['sheep'][ITIME], R)
    for i in np.array(range(Np))[DATA['alive'][ITIME]]:
        dumidx = idx[i][idx[i] != i] # remove identity particle
        if len(dumidx) >= 1: # Some neighbours
            sheepAccTemp = (DATA['sheepVel'][ITIME][DATA['alive'][ITIME]][dumidx].mean(axis = 0) - DATA['sheepVel'][ITIME][i])
            if ITIME < 4.:
                sheepAccFlocking[i] = sheepAccTemp/dt
            else:
                sheepAccFlocking[i] = sheepAccTemp/q*dt
        else:  # No neighbours
            sheepAccFlocking[i] = np.zeros(2)
    return sheepAccFlocking

def doTopo(DATA, ITIME, Np, N):
    tree = KDTree(DATA['sheep'][ITIME][DATA['alive'][ITIME]])
    idx = tree.query(DATA['sheep'][ITIME][DATA['alive'][ITIME]], np.min([N+1, sum(DATA['alive'][ITIME])]))[1][:,1:]
    sheepAccTemp = (DATA['sheepVel'][ITIME][idx].mean(axis = 1) - DATA['sheepVel'][ITIME])
    if ITIME < 4.:
        sheepAccFlocking = sheepAccTemp/dt
    else:
        sheepAccFlocking = sheepAccTemp/q*dt
    return sheepAccFlocking
