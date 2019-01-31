import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint

def f (y, t, params):
	k1, k2 = y
	k, L, xi1, xi2, delta1, delta2, g1, g2, n, T = params
	#print("k1 ",k1, "at t ", t)
	#print("k2 ", k2, " at t ", t)
	derivs = [-(1-xi1*g1(k1,k2,t,n,T))/L + (1-xi2*g2(k1,k2,t,n,T))/L, -(1-xi2*g2(k1,k2,t,n,T))/L + (1-xi1*g1(k1,k2,t,n,T))/L] # out flow - in flow
	return derivs

# Parameters # TODO
k = 40 # average overall density for both rings (should be constant) should be k >= kj/2 to enable multivaluedness
L = 0.25 # length of each ring in miles
#T = 55 # cycle length time (in seconds); # TEST CHANGES IN THIS
T0 = 55
n = 0 # current cycle number
k_j =  60 # jam density ranges from 185-250 per mile per lane typically; they used 60 in example
k_c = k_j / 5 # critical density is approx 1/5 kj
pi1 = 0.45 # green time ratio for ring road 1 # TEST CHANGES IN THIS
pi2 = 0.45 # green time ratio for ring road 2 # TEST CHANGES IN THIS
Delta = (T0 - (T0*(pi1 + pi2)))/2
v_f = 0.01666 # miles per second ; free flow speed 
w_a = 0.01388 # miles per second ; shock-wave speed 

xi1_0 = .85
xi2_0 = .85

def xi1 (t): # the retainment ratio of ring road 1
	val = xi1_0
	return val

def xi2 (t): # the retainment ratio of ring road 2
	val = xi2_0
	return val

def gamma2 (xi1):
	return ((1.0-xi1)*v_f*k_c) / ((L*xi1)*(k_j-k_c))

def delta1(t, n,T): # the signal control based on green ratio, pi1
	signal1 = 0
	if (t >= n*T and t < n*T + pi1*T):
		#print("signal1 is turned on at time ", t)
		signal1 = 1
	return signal1

def delta2(t, n, T): # signal control based on green ratio pi2
	signal2 = 0
	Delta = (T - (T*(pi1 + pi2)))/2
	if(t >= n*T + Delta + pi1*T and t < (n+1)*T - Delta):
		#print("signal2 is turned on at time ",t)
		signal2 = 1

	return signal2

def Q(k_a): # flow-density relationship
	return min(v_f*k_a, w_a*(k_j-k_a))
	
def D(k_a, t): # demand - aka out-ward flow (into junction)
	demand = 0
	if(k_a >= 0 and k_a <= k_c):
		demand = Q(k_a)
	elif (k_a > k_c and k_a <= k_j): # k_a is between k_c and k_j
		demand = Q(k_c)
	else:
		#print("")
		demand = 0
	return demand

def S(k_a, t): # supply - aka in-ward flow (into link from junction)
	supply = 0
	if(k_a >= 0 and k_a <= k_c):
		supply = Q(k_c)
	elif (k_a > k_c and k_a <= k_j): # k_a is between k_c and k_j
		supply = Q(k_a)
	else:
		#print("")
		demand = 0
	return supply

def g1(k1,k2, t,n,T):
	val = delta1 (t,n,T) * min(D(k1,t), S(k1,t)/xi1_0, S(k2,t)/(1-xi1_0))
	#print("g1 at time", t, " is ", val)
	return val

def g2(k1,k2, t,n, T):
	val = delta2 (t,n,T) * min(D(k2,t), S(k2,t)/xi2_0, S(k1,t)/(1-xi2_0))
	#print("g2 at time", t, " is ", val)
	return val

# while loop here which iterates through ranges of values for all the parameters
# initial densities; green ratios; cycle lengths; retainment ratios
# find out patterns, or what are possible values to modify

# for initial settings, determine the region for effective green time for ring 1 -> region 1
# for same initial settings, determine the region for effective green time for ring 2 -> region 2
# each region has defined limits that we can use to determine if the system is in that region or not
# given knowledge of the (k1,k) state, we define (region 1, region 2) such that 
# we can simulate different values of cycle length T, xi1, xi2, and pi1, pi2 to cause the system to transition into the congestion regions 
# the congestion regions are (3,8) and (4,7) in general but this may change according to values of xi1,xi2, and pi1,pi2. It is necessary that 
# xi1, xi2, and pi1, pi2 follow the defined constraints in the paper as well.

# then, we are left with simulating the double ring network model under different parameters and we find
# the ones that cause one or both rings to become congested by diverging it from the original state.

# Since some of these states are stationary and have long-term stability behavior (constant average flow), then it may take a while to determine if the 
# flow truly becomes zero or not.

# We are interested in seeing how changing the signal settings for one cycle can cause negative asymptotic behavior in the traffic flow.

# To determine the flow, we can either map the current average density to the average network flow using the poincare maps or we can simulate for a very long time


for T in range (T0, 500, 10):
	for xi1 in np.arange(xi1_0, 1.0, .1):
		A = (math.exp(gamma2(xi1)*pi1*T) - 1) / (math.exp(gamma2(xi1)*pi1*T) + 1)
		B = (2*k) / (math.exp(gamma2(xi1)*pi1*T)+1)
		print ("T is: ", T, "and xi1 is ", xi1, " and gamma2 is ", gamma2(xi1), " total density 2k is: ", (2*k) )
		print("A is: ", A, " and B is: " , B)
		total = (k_j * A) + B 
		print("total for k1* is: ", total, " and goal is kj: ", k_j)
# 			# Plot the results
# 			fig = plt.figure(1,figsize=(8, 8))
# 			# density vs time
# 			ax1 = fig.add_subplot(311)
# 			ax1.plot(t, soln[:,0])
# 			ax1.set_xlabel('time t')
# 			ax1.set_ylabel('density k1')
# 			# flow vs time
# 			ax2 = fig.add_subplot(312)
# 			ax2.plot(t, soln[:,1])
# 			ax2.set_xlabel('time t')
# 			ax2.set_ylabel('density k2')
# 			# flow vs density - fundamental diagram
# 			ax3 = fig.add_subplot(313)
# 			ax3.plot(soln[:,0], soln[:,1])
# 			ax3.set_xlabel('density k1')
# 			ax3.set_ylabel('density k2')

# for one cycle, change the values of the phase times and determine whether the state is one of the gridlock states

# for init densities
for k1_0 in range(40,60,5):
	# for green ratios
	for k2_0 in range(40, 50, 5):
		# for cycle lengths
		for T in range (T0, 195, 10):
			# Initial Values
			# k1_0 = 40
			# k2_0 = 40
			k = (k1_0 + k2_0)/2
			
			print("k1_0 is: ", k1_0)
			print("k2_0 is: ", k2_0)
			print("k is: ", k)
			print("T is: ", T)
			# Bundle the parameters
			params = [k, L, xi1_0, xi2_0, delta1, delta2, g1, g2, n,T]

			# Bundle the init conditions
			y0 = [k1_0, k2_0]

			# Make time array
			tStop = 150.
			tInc = 0.01
			t = np.arange(0., tStop, tInc)

			# Call the ODE Solver
			soln = odeint(f, y0, t, args=(params,))

			# Plot the results
			fig = plt.figure(1,figsize=(8, 8))
			# density vs time
			ax1 = fig.add_subplot(311)
			ax1.plot(t, soln[:,0])
			ax1.set_xlabel('time t')
			ax1.set_ylabel('density k1')
			# flow vs time
			ax2 = fig.add_subplot(312)
			ax2.plot(t, soln[:,1])
			ax2.set_xlabel('time t')
			ax2.set_ylabel('density k2')
			# flow vs density - fundamental diagram
			ax3 = fig.add_subplot(313)
			ax3.plot(soln[:,0], soln[:,1])
			ax3.set_xlabel('density k1')
			ax3.set_ylabel('density k2')

			plt.tight_layout()
			plt.show()