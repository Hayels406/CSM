import numpy as np

def doFlight(DATA, ITIME, Np, B, ALLSHEEPSEEDOG):
    preyMinusPred = DATA['sheep'][ITIME][DATA['alive'][ITIME]] - DATA['dog'][ITIME]
    normStuff2 = np.linalg.norm(preyMinusPred, axis=1).reshape(sum(DATA['alive'][ITIME]),1)

    Gfunc = B*normStuff2**(-1)

    if ALLSHEEPSEEDOG == 'On':
        DATA['sheepAccFlight'][ITIME][DATA['alive'][ITIME]] = Gfunc*preyMinusPred*(normStuff2**(-1))
    elif ALLSHEEPSEEDOG == 'Off':
        vor = Voronoi(np.concatenate((DATA['sheep'][ITIME][DATA['alive'][ITIME]], DATA['dog'][ITIME][None,:]), axis = 0))
        DATA['dog_index_neighbours'] = np.array([t[1] for t in [(b1, a1) for a1, b1 in vor.ridge_dict.keys()] + vor.ridge_dict.keys() if t[0] == Np]) #Calculates the Voronoi neighbours of dog
        DATA['sheepAccFlight'][ITIME][DATA['alive'][ITIME]][DATA['dog_index_neighbours']] = (Gfunc*preyMinusPred*(normStuff2**(-1)))[DATA['dog_index_neighbours']]
    else:
        sys.exit('Invalid ALLSHEEPSEEDOG:'+ALLSHEEPSEEDOG+' in doAccelerationStep()')
    return DATA
