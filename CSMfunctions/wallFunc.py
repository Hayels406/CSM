import numpy as np

def doWalls(DATA, ITIME, A, B, C, Np, WALLTYPE):
    if WALLTYPE == 'Square':
        for wall in range(len(DATA['walls']['eqn'])):
            w = DATA['walls']['eqn'][wall]

            DATA['forceWalls'][ITIME][DATA['alive'][ITIME]] += (A*np.exp((np.dot(-DATA['walls']['n'][wall], DATA['sheep'][ITIME].T)-np.abs(w[2])+C)/B).reshape(Np,1)*DATA['walls']['n'][wall])[DATA['alive'][ITIME]]

    elif WALLTYPE == 'Circular':
        distance = np.linalg.norm(DATA['sheep'][ITIME][DATA['alive'][ITIME]], axis = 1)
        unit = DATA['sheep'][ITIME][DATA['alive'][ITIME]]/distance.reshape(sum(DATA['alive'][ITIME]),1)
        DATA['forceWalls'][ITIME][DATA['alive'][ITIME]] = -A*unit*np.exp((distance-DATA['walls'][0]+C)/B ).reshape(sum(DATA['alive'][ITIME]),1)
    return DATA
