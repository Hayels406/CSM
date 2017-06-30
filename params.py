TF = 2000.
dt = 0.1
NP = 200
a = 0.01
b = 0.1
p = 3.0
c = 1.5
f = 1.
A = 2*10**3
B = 0.08
epsilon = dt**6
sheepSize = 0.45
sheepMass = 60.
dogMass = 17.

r = 10.
n = 100
visualDist = 20.
tau = 1.
accCap = 1.
noSegments = 6

timeStepMethod = 'Adaptive'   #Euler or Adaptive
updateMethod = 'Acceleration' #Acceleration or Velocity
flocking = 'PredatorPrey'           #PredatorPrey, Vicsek, TopoVicsek

sheepSheepInteration = 'On' #F_ij term
allSheepSeeDog = 'On' #Can all sheep see dog?
gaussian = 'Off'
normDog = 'Off'
cap = 'Off'
segments = 'On'

wallType = 'Square'
wallLeft = -50
wallRight = 50
wallBottom = -50
wallTop = 50

snapshotPeriod = 5000.0 #units of time
plotPeriod = 50. #units of time
plot = 'On' #On or Off
savePlotPng = 'Off' #On or Off
saveDataH5 = 'On'

dog_init = [-5., -10.]		    #np.random.rand(2)*3.0+1.
sheep_init = 'Random'		#Random or Grid
dog_vel_init = [0., 0.]
sheep_vel_init = 0.1
sheep_area = 80.

loadFromFile = 'Off'
fileName = 'data-0000200.h5'

#plotting tweaks
showSheepDogCanSee = 'Off'
showDogInfluence = 'Off'
segmentColours = 'On'
