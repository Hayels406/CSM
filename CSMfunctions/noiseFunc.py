import numpy as np
import sys

def includeNoise(DATA, ITIME, NOISE):
    if (NOISE == 'Uniform')*(ITIME > 4):
        dogTheta = np.arctan2(DATA['dogVel'][ITIME,1], DATA['dogVel'][ITIME,0])
        sheepTheta = np.arctan2(DATA['sheepVel'][ITIME,:,1][DATA['alive'][ITIME]], DATA['sheepVel'][ITIME,:,0][DATA['alive'][ITIME]])
        dogVel = np.sqrt((DATA['dogVel'][ITIME]**2).sum())
        sheepVel = np.sqrt((DATA['sheepVel'][ITIME][DATA['alive'][ITIME]]**2).sum(axis = 1))
        dogTheta = dogTheta + np.random.rand(1)[0]*np.pi*eta - np.pi*eta/2.
        sheepTheta = sheepTheta + np.random.rand(sum(DATA['alive'][ITIME]))*np.pi*eta - np.pi*eta/2.
        DATA['dogVel'][ITIME] = dogVel*np.array([np.cos(dogTheta), np.sin(dogTheta)])
        DATA['sheepVel'][ITIME][DATA['alive'][ITIME]] = sheepVel.reshape(sum(DATA['alive'][ITIME]), 1)*np.transpose(np.array([np.cos(sheepTheta), np.sin(sheepTheta)]))
    elif (NOISE == 'Normal')*(ITIME > 4):
        dogTheta = np.arctan2(DATA['dogVel'][ITIME,1], DATA['dogVel'][ITIME,0])
        sheepTheta = np.arctan2(DATA['sheepVel'][ITIME,:,1][DATA['alive'][ITIME]], DATA['sheepVel'][ITIME,:,0][DATA['alive'][ITIME]])
        dogVel = np.sqrt((DATA['dogVel'][ITIME]**2).sum())
        sheepVel = np.sqrt((DATA['sheepVel'][ITIME][DATA['alive'][ITIME]]**2).sum(axis = 1))
        dogTheta = dogTheta + np.random.normal(0, sigma)
        sheepTheta = sheepTheta + np.random.normal(0, sigma, sum(DATA['alive'][ITIME]))
        DATA['dogVel'][ITIME] = dogVel*np.array([np.cos(dogTheta), np.sin(dogTheta)])
        DATA['sheepVel'][ITIME][DATA['alive'][ITIME]] = sheepVel.reshape(sum(DATA['alive'][ITIME]), 1)*np.transpose(np.array([np.cos(sheepTheta), np.sin(sheepTheta)]))
    elif NOISE != 'Off':
        sys.exit('Invalid Noise')
    return DATA
