###########################################################################################################################################
### Shklarsh Model ### 
#Sklarsh model proposed heuristics for simulating bacterial food search using a gradient decent algorithm. The assumptions are as follows:
#1. Individually Bacteria move towards higher food concentration with a small randomness attached, which depends on difference in food concentration(gradient) 
#2. Collectively there are three interactions:
#2.1 Repulsion - cells too close are repulsed to avoid collision
#2.2 Orientaion - to match direction of neighbouring cells
#2.3 Attraction - from distant agents to keep the group unified
#3. Cells can identify which cells lie in which range and update thier velocity vectors accordingly 

###########################################################################################################################################
### Distributed Gradient Decent ###
#DGD is a type of gradient decent which uses hueristics to closly simulate bacterial food searches. Many changes were proposed in this model to make it 
#more "realistic" than its predecessor, the Sklarsh model. The difference lies in the cell-cell interaction. The model assumes that cells cannot identify precisly which 
#radius do other cells lie in and the whole network of cells is not static. Interactions are divided into physical (repulsion) and chemical (orientation, attraction).
#The individual movement and repulsion interations are same as Sklarsh model. Orientation and attraction coefficents are introduced which dictate the influence of all the cells on
#on the i-th cell. This system relies on very small message sizes over the network of the sizes of (logL + 3) bits. Where L is the number of discrete values that the magnitude
#of velocity vector can take and 3 comes from log8, where 8 is the number of cardinal directions. 
#DGD can be considered to be a more generalised form of Sklarsh model where the attraction and orientation radius are 
#not quantised but diffussed over the whole terrain.
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.spatial import distance

# global parameters
pi = np.pi
food = np.matrix([3, 2]) # global maxima
agents = 30
alpha = .001
kk = 35

def delta_t(d):
    return np.exp(alpha*(np.square(d)-np.square(kk)))    # variable step size
#delta_t = 0.5                                             # constant step size
num_of_iterations = 400
initial_d = 20
found_radius = 2.5
sigma = 0.01
K = 200
variance = 1000
RR = 0.1 # repulsion radius
RO = 4   # oreintaion radius
RA = 4.3 # attraction radius
L = 4    # descrete levels
T = 1    # threshold for descretization function
collision_delay = 0
#-- attraction and orientation functions, can be changed to zero for no communication
def orientFunc(x):
    x = np.exp(-2*x)
    return x
def attractFunc(x):
    x = np.exp(-0.5*x)
    return x

#-- Objective function (terrain)
def nf(x,y):
	return math.sqrt(x**2 + y**2)

def obstacles(x,y): # local maxima
	return min(0, -4*math.cos(pi*x/3)*math.cos(pi*y/3) + 0.5)

def foodSource(x,y): # gradient function
	return (K / (np.sqrt(2*pi*variance))) * np.exp(-(nf(x - food[0,0] ,y - food[0,1])**2 / (2*variance) ))

def terrain(x,y):
	return foodSource(x,y) + obstacles(x,y)

#-- The swarm matrix(somewhere on the terrain)
theta = np.random.random()*pi*2
x_min = math.floor(initial_d*math.cos(theta)) + food[0,0] - 0.5
y_min = math.floor(initial_d*math.sin(theta)) + food[0,1] - 0.5
X_new = 3*np.random.random_sample([agents, 1]) + x_min
Y_new = 3*np.random.random_sample([agents, 1]) + y_min
X = np.concatenate((X_new, Y_new), axis = 1)# (agents X 2) matrix 


Y = np.concatenate((X,food))
plt.scatter(Y[:,0], Y[:,1])
plt.scatter(food[0,0], food[0,1], s = pi*(5**2), color = 'g')
plt.axis([-25,25,-25,25])
plt.savefig("before.png")
plt.close()


# initialize velocity matrix
V = np.zeros((agents, 2))

#initialize Collision delay matrix
KO = np.zeros((agents,1))

# initialize Timeline matrix
converged_agents_timeline = np.zeros(num_of_iterations/10)

# weigth function
def wf(u,v, delta_gradient):
    return u + 30*v

# check quadrant for u vector and modify value of theta accordingly
def check_quadrant_matrix(V_i, theta):
    x = V_i[0,0]
    y = V_i[0,1]
    if x > 0 and y < 0:
	return theta + pi
    elif x < 0 and y < 0:
	return theta + pi
    return theta

# descritization function, this implementation is based on the theory in the paper, i could not find any mathematical formula for exact numbers. 
def df(x, T=T):
    # calculate magnitude of u
    abs_x = np.sqrt(x[0,0]**2 + x[0,1]**2)
    # calculate direction of u
    if abs_x > 0:
	unit_x = x/abs_x
    else:
	unit_x = x   

    # for quantizing direction of u
    direction = [0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4, 3*pi/2, 7*pi/4] 
    limits = [pi/8, 3*pi/8, 5*pi/8, 7*pi/8, pi+(pi/8), pi+(3*pi/8), pi+(5*pi/8), pi+(7*pi/8)]

    if abs_x == 0:
	hypotenous = 1
    else:
	hypotenous = abs_x
    theta = np.arccos(x[0,0]/hypotenous)
    theta = check_quadrant_matrix(x, theta)

    if limits[7] < theta <= limits[0]:
	theta = direction[0]
    for i in range(7):
        if limits[i] < theta <= limits[i+1]:
	    theta = direction[i+1]
	    break
    unit_x = np.matrix(([np.cos(theta), np.sin(theta)]))

    # for quantizing direction of u
    level = [0.25, 0.50, 0.75, T]

    if abs_x >=T: 
        abs_x = level[3]
    elif level[2] <= abs_x < level[3]:
	abs_x = level[2]
    elif level[1] <= abs_x < level[2]:
	abs_x = level[1]
    elif level[0] <= abs_x < level[1]:
	abs_x = level[0]
    else:
	abs_x = 0
    
    return abs_x*unit_x 
    

# check quadrant for velocity vector and modify value of theta accordingly
def check_quadrant(V_i, theta):
    x = V_i[0]
    y = V_i[1]
    if x > 0 and y < 0:
	return theta + pi
    elif x < 0 and y < 0:
	return theta + pi
    return theta

# distance function to calculate number of agents that have found food source
def distFunc(X, food, found_radius, i):
    distanceMatrix = np.sqrt(np.sum(np.square(np.subtract(X,food)),axis = 1))
    flag_stop = False  
    for j in range(np.size(distanceMatrix)):
	if distanceMatrix[j,0] < found_radius:
	    distanceMatrix[j,0] = 1
	else:
	    distanceMatrix[j,0] = 0
    converged_agents = np.sum(distanceMatrix)
    if i%10 == 0:
        converged_agents_timeline[i/10] = converged_agents 
    if converged_agents >= (3*agents/4):
        flag_stop = True
    return converged_agents_timeline, flag_stop

#-- interaction function to account for effect of other agents on i-th agent
def interaction(i, X, V, D_i):
    u = np.zeros((1,2))
    flag_RR = 1
    # check if any agent is in Repulsion Radius, if yes, take their sum and calculate u
    for j in range(agents):
        if D_i[j] <= RR and j != i:
	    u += -(X[j,:] - X[i,:])
	    flag_RR = 0
    if np.linalg.norm(u) > 0:
	u = u/np.linalg.norm(u)
    # if no agents in repulsion radius then calculate total orientation and attraction and return u
    if flag_RR:
#	DistanceArray = np.sqrt(np.sum(np.square(np.subtract(X,X[i,:])),axis = 1))
#	w_o = np.tile(np.reshape(orientFunc(DistanceArray),(agents,1)),2)
#	w_a = np.tile(np.reshape(attractFunc(DistanceArray),(agents,1)),2)
	w_o = np.tile(orientFunc(D_i),2)
	w_a = np.tile(attractFunc(D_i),2)
	w_o[i] = 0 # no interaction with itself
	w_a[i] = 0 # no interaction with itself

	attracts = X - X[i,:]
	
	orient = np.zeros((1,2))
	attract = np.zeros((1,2))
	for j in range(agents):
	    orient += np.multiply(V[j,:], w_o[j])
	    abs_attracts = np.sqrt(attracts[j,0]**2 + attracts[j,1]**2)
	    if abs_attracts == 0:
	        attracts[j,:] = 0 
	        attract += np.multiply(attracts[j,:], w_a[j])
	    else:
		attract += np.multiply(attracts[j,:]/abs_attracts, w_a[j])	        

	u = df(orient + attract)

    return u  
	
    

#-- velocity function for updation V_new
def velocity(i, X, V, D_i):
    DistanceArray = np.sqrt(np.sum(np.square(np.subtract(X[i,:],food)),axis = 1))  # comment out if using constant time step
    prevX = X[i,:] - delta_t(DistanceArray)*V[i,:]                                 # comment out if using constant time step 
#    prevX = X[i,:] - delta_t*V[i,:]  # use for constant time step
    delta_gradient = terrain(X[i,0],X[i,1]) - terrain(prevX[0,0],prevX[0,1])       # delete axis 0 from prevX terms if using constant time step
    hypotenius = np.sqrt(np.square(V[i,0]) + np.square(V[i,1]))
    if hypotenius == 0:
	hypotenius = 1
    theta = np.arccos(V[i,0]/hypotenius) 
    mu = check_quadrant(V[i,:], theta)
    if delta_gradient <= 0:
        theta_new = np.random.normal(mu, pi)
    else:
	theta_new = mu
    V[i,:] = [np.cos(theta_new), np.sin(theta_new)] 
    v = wf(interaction(i, X, V, D_i), V[i,:], delta_gradient)
# normalize the value of v
    if np.linalg.norm(v) > 0:
	v = v/np.linalg.norm(v)
# add collion delay for agents with pairwise distances < RR
    if (KO[i,0] > 0):
        v = 0
	KO[i,0] -= 1
    for j in range(np.size(D_i)):
	if D_i[j] > 0 and D_i[j] < RR:
	    KO[i,0] = collision_delay
	    break
# add noise	    
    v = v + np.random.random_sample([1,2])*sigma
    return v

#-- velocities function for updating V
def velocities(V,X):
    V_new = np.zeros((agents,2))
# for pair wise distance between agents to calculate interations later on
    D = distance.cdist(X, X, 'euclidean')
# calculate velocity for each agent in matrix X
    for i in range(agents):
	V_new[i,:] = velocity(i, X, V, D[:,i])
    return V_new

#-- start simulation
for i in range(num_of_iterations):
    V = velocities(V,X)
    DistanceArray = np.sqrt(np.sum(np.square(np.subtract(X,food)),axis = 1)) # comment out if using constant time step
    X = X + np.multiply(delta_t(DistanceArray),V)                           # comment out if using constant time step
#    X = X + delta_t*V                                                       # use for constant step size
    converged_agents_timeline, flag_stop = distFunc(X, food, found_radius,i)
    if flag_stop:
	print "Agents have reached food source after %d iterations" % i
	break

    	
print converged_agents_timeline

#-- plot the agents and food
plt.scatter(food[0,0], food[0,1], s = pi*(5**2), color = 'g')
plt.scatter(X[:,0], X[:,1])	
plt.axis([-25,25,-25,25])
plt.savefig("after.png")
plt.close()
#-- plot slope of model
iteration_axis = np.linspace(0,400,40)
plt.plot(iteration_axis, converged_agents_timeline)
plt.xlabel("iterations")
plt.ylabel("No of agents converged")
plt.axis([0 ,400,0 ,30])
plt.savefig("timeline.png")
plt.close()
