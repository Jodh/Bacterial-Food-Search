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

# global parameters
pi = np.pi
food = np.matrix([3, 2]) # global maxima
agents = 30 
delta_t = 0.5 
num_of_iterations = 100
initial_d = 20

#-- Objective function (terrain)
def nf(x,y):
	return math.sqrt(x**2 + y**2)

def obstacles(x,y): # local maxima
	return min(0, -4*math.cos(pi*x/3)*math.cos(pi*y/3) + 0.5)

def foodSource(x,y): # gradient function
	return min(2.5, -2*log(nf(x-3, y-2)**2))

def terrain(x,y):
	return foodSource(x,y) + obstacles(x,y)

#-- The swarm matrix(somewhere on the terrain)
theta = np.random.random()*pi*2
x_min = floor(initial_d*math.cos(theta)) + food[0,0] - 0.5
y_min = floor(initial_d*math.sin(theta)) + food[0,1] - 0.5
X_new = 3*np.random.random_sample([agents, 1]) + x_min
Y_new = 3*np.random.random_sample([agents, 1]) + y_min
X = np.concatenate((X_new, Y_new), axis = 1)# (agents X 2) matrix 
print X

#-- update function for agent positions
def update:
	V = zero[size(X)] #initialize velocity matrix 
	KO = zero[agent_number,1] #initialize collision_delay_matrix 
	for i in range(interation):
		if half of the agents reach [c1, c2]: break
		#-- update V matrix
		D = pair_wise_distance_matrix #Distance matrix D for every pair of agents. n X n matrix
		
		for i in range(1,(agents +1)):
			# calculate gradient change for i-th agent, delta_g 
			X_prev(i) = X(i) - dt*V(i) 	
			delta_g = terrain(X(i)) - terrain(X_prev(i))
			# tumble if delta_g <= 0
			if delta_g <= 0: pick random travel direction
			else keep going in same direction
			#-- calculate orentaion and attraction components for i-th agent over all agents
			# |-equation (7) from paper	
			#-- do not calculate V for agents where collision is detected by checking Repulsion Radius and D(i) matrix 
		#-- Update X matrix 
		X = X + delta_t*V 	
