import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
from scipy.integrate import odeint
import sys 

pi1_init = float(input('What is pi1_init?   '))
pi2_init = float(input('What is pi2_init?   '))
xi1 = float(input('What is xi1?   '))
xi2 = float(input('What is xi2?   '))


print(pi1_init)
print(pi2_init)
print(xi1)
print(xi2)
# pi1_init = 0.5
# pi2_init = 0.5

pi1_delta = 0.1
pi2_delta = 0.1
T = 55 # seconds
#xi1 = .55
#xi2 = .55
k_j = 150 # jam density ranges from 185-250 per mile per lane typically; they used 60 in example
k_c = k_j / 5 # critical density is approx 1/5 kj
k_1_init = 50
k_2_init = 100
k = 75 # average overall density for both rings (should be constant) should be k >= kj/2 to enable multivaluedness
v_f = 0.01666 # miles per second ; free flow speed 
L = 0.25 # length of each ring in miles
# constraints
#pi1 + pi2 = 1
#xi1, xi2 > 0.5
pi1_values = []
pi2_values =[]
k1_values = set([])
xi1_values = []
xi2_values = []

xi1_values = np.arange(0.51, .99, 0.02)
xi2_values = np.arange(0.51, .99, 0.02)
pi1_deltas = np.arange(pi1_init-1, 1-pi1_init , 0.01)

def dk1_dt_r2(k_1, t): 
	return (t * -1) * (1-xi1) * v_f * k_c / L 
def dk1_dt_r3 (k_1, t):
	return (t * -1 * (k_j-k_1) * ((1-xi1)*v_f*k_c) / (L*xi1*(k_j-k_1)) ) 
def dk1_dt_r6 (k_1, t):
	return (t * (1-xi2) * v_f * k_c / L ) 


def region_2_check(k_1,k):
	in_region = 0
	if ((k_1/2 <= k) and (k <= (xi1*k_j + (1-xi1)*k_c + k_1)/2) and (k_c <= k_1) and (k_1 < k_j-xi1*(k_j-k_c))):
		in_region = 1
	else:
		in_region = 0
	return (in_region)
#Aho
def region_3_check (k_1,k):
	in_region = 0
	if ((k_1/2 <= k) and (k <= (2*xi1-1)/(2*xi1)*k_j + k_1/(2*xi1)) and ((k_j - xi1* (k_j-k_c)) <= k_1) and (k_1 <= k_j)):
		in_region = 1
	else:
		in_region = 0
	return (in_region)

def region_6_check(k_1,k):
	in_region = 0
	if (((k_1+k_c)/2 < k) and ((k_j+k_1-(xi2*(k_j-k_c)))/2) and k_1 >= 0 and (k_1 <= k_j - (1-xi2)*(k_j-k_c))):
		in_region = 1
	else:
		in_region = 0
	return (in_region)

def region_8_check (k_1,k) :
	in_region = 0
	if ((k_1/2 + (k_c*(k_j-k_1))/(2*(1-xi2)*(k_j-k_c)) < k) and (k < ((1-2*xi2)*k_j + k_1)/(2*(1-xi2))) and ((k_j-(1-xi2)*(k_j-k_c)) < k_1) and (k_1 <= k_j) ):
		in_region = 1
	else:
		in_region = 0
	return (in_region)

# for i in range(1, 20):
# 	k_1_init  = k_1_init + 5 
# 	k_2_init = k * 2 - k_1_init
# 	if (region_3_check(k_1_init) == 0):
# 		print("k_1 is NOT in region 3")
# 		print(k_1_init)
# 	else:
# 		print("k_1 is IN region 3")
# 		print(k_1_init)

# we will test across various values for k1,k
k1min = round(k_j/2)
k1max = round(k_j)
# k must be greater than region 8 max and lower than region 6 minimum i.e. region 3 ^ region 8
kmin = round(k_j/2)
kmax = round((xi1*k_j+(1-xi1)*k_c+k1max) /2)

# for k1 in range (k1min,k1max,1):
# 	for k in range (kmin, kmax, 1):
# 		if(region_3_check(k1,k)):
# 			print("In region 3 ", k1, " ", k)
# 		else: 
# 			print("Not in region 3 ", k1, " ", k)
# xi1 = .55
# xi2 = .95
# print(k1min)
# print(k1max)
# print(kmin)
# print(kmax)
pi1_deltas2 = []
xi1_values_valid = []
xi2_values_valid = []
xi1_values2 = set([])
xi2_values2 = set([])
# first check to see which k1,k points satisfy 2,6
for k1 in range(k1min, k1max, 1):
	for k in range (kmin, kmax, 1):
		# for xi1 in xi1_values:
		# 	for xi2 in xi2_values:
		for pi1_delta in pi1_deltas:
			pi1 = pi1_init + pi1_delta
			pi2 = 1 - pi1
			#k1 must be in region 2 at the start of phase 1 of first cycle
			#print("Before region2: ",k1)
			if(region_2_check(k1,k) == 0):
				#print("Not in region2")
				continue
			k1_new = max(0, k1 + dk1_dt_r2(k1, pi1*T))
			#print("After region2: ",k1_new)
			#k1 must in region 6 at the start of phase 2 of first cycle
			if(region_6_check(k1_new,k) == 0):
				#print("Not in region6")
				continue
			k1_new = k1_new + dk1_dt_r6(k1_new, pi2*T)
			#print("After region6: ",k1_new)
			#k1 must be in region 2 at the start of phase 1 of next cycle
			if(region_2_check(k1_new,k) == 0):
				#print("Not in region 2 in following cycle")
				continue
			# if a point satisfies these conditions and if k is in region 3 as well, we should store the point
			if ((k1_new/2 <= k) and (k <= (2*xi1-1)*k_j/(2*xi1)+ k1_new/(2*xi1))):
				k1_values.add(round(k1_new))
				xi1_values2.add(xi1)
				xi2_values2.add(xi2)
				pi1_deltas2.append(pi1_delta)

k_values = np.arange(kmin, kmax, 1)
k1_values_valid = list(set(k1_values))
xi1_values_valid = list(set(xi1_values2))
xi2_values_valid = list(set(xi2_values2))
pi1_deltas_valid = list(set(pi1_deltas2))
print(k_values)
print(k1_values_valid)
print(xi1_values_valid)
print(xi2_values_valid)
print(pi1_deltas_valid)
if(k_values==[] or k1_values_valid==[]):
	print("Empty array")
	exit(0)

x=[]
y=[]
z=[]
xi1s = []
xi2s = []
# xi1 = .55
# xi2 = .50
# then check whether k1,k that are in 2,6 can go to region 3 or 8 to make congestion after changing pi values
pi1_min = 0
pi1_d_min = 0
z2 = []
k1_finals = set([])
for k_1 in k1_values_valid:
	for k in k_values:
		#for xi1 in xi1_values_valid:
		#	for xi2 in xi2_values_valid:
		updated = 0
		k1_min = 10000
		pi1_min = -10000
		pi1_d_min = 100000 
		# determine which pi1 value for the current k,k1 state 
		# minimizes the effect on the k,k1 state yet causes successful attack
		for pi1_delta in pi1_deltas:
			pi1 = pi1_init + pi1_delta
			pi2 = 1.0 - pi1
			#k1 after a cycle; should be in region 3
			#print("xi1, xi2: ", xi2, " ", xi2)
			#print("pi1, pi2: ", pi1, " ", pi2)
			# print("Before phase 1: " , k_1, " ", k, " ", pi1)
			# print("dk1dtr2: ",dk1_dt_r2(k_1, pi1*T))
			k1_new = k_1 + round(max(0,dk1_dt_r2(k_1, pi1*T)))
			# print("After phase 1: ", k1_new, " ", k, " ", pi1)
			k1_new = k1_new + round(dk1_dt_r6(k1_new, pi2*T)) 
			#print("After 2 phases: ", k1_new, " ", k, " ", pi1)
			# check if in region 3 now
			#if(k1_new > k_j - (1-xi2)*(k_j-k_c)):
			if(region_3_check(k1_new,k)==1):
				# add this to the x plot points
				#if(k1_new < k1_min):
				if(abs(pi1_delta) < abs(pi1_d_min) and pi1 > 0.2 and pi1 < 0.8):
					k1_min = k_1
					pi1_min = pi1
					pi1_d_min = pi1_delta
					updated = 1
					k1_final = k1_new
					#print("In region 3", k1_new, " ", k)
			# else:
			# 	print("Not in region 3", k1_new, " ", k)

		# add the k,k1,pi1 values
		if(updated == 1):
			k1_finals.add(k1_new)
			x.append(k1_min)
			y.append(k)
			z.append(pi1_d_min) # contains all the optimal pi1 values for each k,k1
			z2.append(pi1_min)
			#xi1s.append(xi1)
			#xi2s.append(xi2)
final_x = set(x)
final_y = set(y)
final_z = set(z)
final_z2 = set(z2)
print(final_x)
print(k1_finals)
print(final_y)
print(final_z)
print(final_z2)
#print(xi1s)
#print(xi2s)

# # Plot the results
fig = plt.figure(1,figsize=(8, 8))
# # density vs time
ax1 = fig.add_subplot(211, projection='3d')
ax1.plot(x, y, z)
ax1.set_xlabel('k1')
ax1.set_ylabel('k')
ax1.set_zlabel('pi1_delta')
ax1.set_title(['pi1: %f, pi2: %f, xi1: %f, and xi2: %f' %(pi1_init, pi2_init, xi1, xi2)])
ax2 = fig.add_subplot(212, projection='3d')
ax2.plot(x,y,z2)
ax2.set_xlabel('k1')
ax2.set_ylabel('k')
ax2.set_zlabel('pi1_modified')
ax2.set_title(['pi1: %f, pi2: %f, xi1: %f, and xi2: %f ' %(pi1_init, pi2_init, xi1, xi2)])
plt.show()
