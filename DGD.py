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
delta_t = 0.5 
num_of_iterations = 100
initial_d = 20
found_radius = 0.9
sigma = 0.1
RR = 0.1 # repulsion radius
RO = 4   # oreintaion radius
RA = 4.3 # attraction radius
L = 4    # descrete levels
T = 1    # threshold for descretization function
collision_delay = 10
#-- attraction and orientation functions, can be changed to zero for no communication
def orientFunc(x):
    x = x*1
    return x
def attractFunc(x):
    x = x*1
    return x

#-- Objective function (terrain)
def nf(x,y):
	return math.sqrt(x**2 + y**2)

def obstacles(x,y): # local maxima
	return min(0, -4*math.cos(pi*x/3)*math.cos(pi*y/3) + 0.5)

def foodSource(x,y): # gradient function
	return min(2.5, -2*math.log(nf(x-3, y-2)**2))

def terrain(x,y):
	return foodSource(x,y) + obstacles(x,y)

#-- The swarm matrix(somewhere on the terrain)
theta = np.random.random()*pi*2
x_min = math.floor(initial_d*math.cos(theta)) + food[0,0] - 0.5
y_min = math.floor(initial_d*math.sin(theta)) + food[0,1] - 0.5
X_new = 3*np.random.random_sample([agents, 1]) + x_min
Y_new = 3*np.random.random_sample([agents, 1]) + y_min
X = np.concatenate((X_new, Y_new), axis = 1)# (agents X 2) matrix 

#x_food = np.matrix([food[0,0]])
#y_food = np.matrix([food[0,1]])
#X_new = np.concatenate((X_new,x_food), axis = 0)  
#Y_new = np.concatenate((Y_new,y_food), axis = 0)
plt.scatter(X_new, Y_new)
plt.show()


# initialize velocity matrix
V = np.zeros((agents, 2))

#initialize Collision delay matrix
KO = np.zeros((agents,1))


# weigth function
def wf(u,v,delta_gradient):
    return u + 10*(delta_gradient > 0)*v

# descritization function
def df(x):
    return x

# distance function to calculate number of agents that have found food source
def distFunc(X, food, found_radius):
    distanceMatrix = np.sqrt(np.sum(np.square(np.subtract(X,food)),axis = 1))
    for i in range(np.size(distanceMatrix)):
	if distanceMatrix[i,0] < found_radius:
	    distanceMatrix[i,0] = 1
	else:
	    distanceMatrix[i,0] = 0
    if np.sum(distanceMatrix) >= (agents/2):
	return True
    else:
	return False

#-- interaction function to account for effect of other agents on i-th agent
def interaction(i, X, V, D_i):
    repulsion_index = np.where((D_i <= RR) & (D_i > 0))
    if collision_delay == 0 and np.any(repulsion_index):
	u = -1*np.sum(X[repulsion_index,:] - X[i,:], axis = 0)
    else:
	w_o = np.tile(np.reshape(orientFunc(D_i),(agents,1)),2)
	w_a = np.tile(np.reshape(attractFunc(D_i),(agents,1)),2)
	w_o[i] = 0 # no interaction with itself
	w_a[i] = 0 # no interaction with itself

	attracts = X - X[i,:]
	
	orient = np.zeros((1,2))
	attract = np.zeros((1,2))
	for j in range(agents):
	    orient += df(V[j,:]*w_o[j])
	    attract += df(attracts[j,:]*w_a[j])

	u = orient + attract
	if np.linalg.norm(u) > 0:
	    u = u/np.linalg.norm(u)
    return u  
	
    

#-- velocity function for updation V_new
def velocity(i, X, V, D_i):
    prevX = X[i,:] - delta_t*V[i,:]
    delta_gradient = terrain(X[i,0],X[i,1]) - terrain(prevX[0],prevX[1])
    theta = pi*np.random.random()*(delta_gradient <=0)
    R = np.matrix([[np.cos(theta), -1*np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    v = wf(interaction(i, X, V, D_i), V[i,:]*R, delta_gradient)
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


for i in range(num_of_iterations):
    V = velocities(V,X)
    X = X + delta_t*V
    if distFunc(X, food, found_radius):
	print X, i, "agents have reached food source"
	break

distanceMatrix = np.sqrt(np.sum(np.square(np.subtract(X,food)),axis = 1))
print distanceMatrix
plt.scatter(X[:,0], X[:,1])
plt.show()

    
