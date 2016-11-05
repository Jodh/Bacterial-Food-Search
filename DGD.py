###########################################################################################################################################
### Shklarsh Model ### 
"""Sklarsh model proposed heuristics for simulating bacterial food search using a gradient decent algorithm. The assumptions are as follows:
1. Individually Bacteria move towards higher food concentration with a small randomness attached, which depends on difference in food concentration(gradient) 
2. Collectively there are three interactions:
2.1 Repulsion - cells too close are repulsed to avoid collision
2.2 Orientaion - to match direction of neighbouring cells
2.3 Attraction - from distant agents to keep the group unified
3. Cells can identify which cells lie in which range and update thier velocity vectors accordingly""" 

###########################################################################################################################################
### Distributed Gradient Decent ###
""" DGD is a type of gradient decent which uses hueristics to closly simulate bacterial food searches. Many changes were proposed in this model to make it 
more "realistic" than its predecessor, the Sklarsh model. The difference lies in the cell-cell interaction. The model assumes that cells cannot identify precisly which 
radius do other cells lie in and the whole network of cells is not static. Interactions are divided into physical (repulsion) and chemical (orientation, attraction).
The individual movement and repulsion interations are same as Sklarsh model. Orientation and attraction coefficents are introduced which dictate the influence of all the cells on
on the i-th cell. This system relies on very small message sizes over the network of the sizes of (logL + 3) bits. Where L is the number of discrete values that the magnitude
of velocity vector can take and 3 comes from log8, where 8 is the number of cardinal directions. 
DGD can be considered to be a more generalised form of Sklarsh model where the attraction and orientation radius are 
not quantised but diffussed over the whole terrain."""

#-- Objective function (terrain)
food = guassain distribution, with center at [c1, c2]
obstacles = half cosine grid  
terrain = food + obstacles
 
#-- Upper and lower bounds for swarm
-20 <= x_values, y_values <= 20

#-- The swarm matrix(somewhere on the terrain)
agent_number = n 
X = randomly assign angle and radian distance from [c1, c2] # (agent_number X 2) matrix 

#-- update function for agent positions
def update:
	V = zero[size(X)] #initialize velocity matrix 
	KO = zero[agent_number,1] #initialize collision_delay_matrix 
	for i in range(interation):
		if half of the agents reach [c1, c2]: break
		#-- update V matrix
		D = pair_wise_distance_matrix #Distance matrix D for every pair of agents. n X n matrix
		
		for i in range(1,(agent_number +1)):
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
		X = X + dt*V 	
