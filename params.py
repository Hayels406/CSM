TF = 100.
dt = 0.1
NP = 400
a = 1.
b = 0.1
p = 3.0
c = 10.
f = 1.

sheepMass = 0.5
dogMass = 1.

r = 10.
groupSize = 2
n = 50
visualDist = 20.
tau = 1.
accCap = 1.
noSegments = 6

timeStepMethod = 'Adaptive'         #Euler or Adaptive
updateMethod = 'Acceleration'       #Acceleration or Velocity
flocking = 'PredatorPrey'           #PredatorPrey, Vicsek, TopoVicsek

sheepSheepInteration = 'Off'    #F_ij term
allSheepSeeDog = 'On'           #Can all sheep see dog?
gaussian = 'Off'
normDog = 'Off'                 #Normalise the dogs acceleration
cap = 'Off'                     #Cap the acceleration of the sheep
segments = 'Off'                #Dog only interacts with sheep in the segment in front of it
emergencyCheck = 'Off'

wallType = 'Square'           #Square or Circular
wallLeft = -1.25
wallRight = 1.25
wallBottom = -1.25
wallTop = 1.25

snapshotPeriod = 100.0 #units of time
plotPeriod = 0.05 #units of time
plot = 'On' #On or Off
savePlotPng = 'Off' #On or Off
saveDataH5 = 'Off'

dog_init = [-0., -0.]		    #np.random.rand(2)*3.0+1.
sheep_init = 'Random'		#Random or Grid
dog_vel_init = [0., 0.]
sheep_vel_init = 0.
sheep_area = 1.

loadFromFile = 'Off'
fileName = 'data-0000200.h5'

#plotting tweaks
showSheepDogCanSee = 'Off'
showDogInfluence = 'Off'
segmentColours = 'Off'
