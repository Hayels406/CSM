TF = 5000.
dt = 0.1
NP = 200
a = 0.01
b = 1.0
p = 3.0
c = 1.5
f = 1.
A = 2*10**3
B = 0.08
epsilon = dt**6
sheepSize = 0.45
sheepMass = 60.

r = 10.
n = 199
visualDist = 20.
tau = 1.
accCap = 1.
noSegments = 20

timeStepMethod = 'Adaptive'   #Euler or Adaptive
updateMethod = 'Acceleration' #Acceleration or Velocity
flocking = 'PredatorPrey'           #PredatorPrey, Vicsek, TopoVicsek

sheepSheepInteration = 'On' #F_ij term
allSheepSeeDog = 'On'       #Can all the sheep see the dog
gaussian = 'Off'            #New f function
normDog = 'Off'
cap = 'Off'
segments = 'Off'

wallType = 'Square'
wallLeft = -50
wallRight = 50
wallBottom = -50
wallTop = 50

snapshotPeriod = 1000.0  #units of time
plotPeriod = 50.        #units of time
plot = 'Off'            #On or Off
savePlotPng = 'On'      #On or Off
saveDataH5 = 'On'       #On or Off

dog_init = [-5., -10.]		#np.random.rand(2)*3.0+1.
sheep_init = "Random"     #Grid or Random
dog_vel_init = [0., 0.]   #Initial velocity of the dog
sheep_vel_init = 0.1      #Initial speed of sheep when Random
sheep_area = 80.          #Initial size of sheep when Random

loadFromFile = 'Off'
fileName = 'data-0001000.h5'

#Plotting tweaks
showSheepDogCanSee = 'Off'
showDogInfluence = 'Off'
segmentColours = 'Off'
