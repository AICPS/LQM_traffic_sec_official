import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
from scipy.integrate import odeint
import sys 
from interval import interval, inf, imath, fpu

def graph(formula, x_range, k_1):  
    x = np.array(x_range)  
    y = formula(k_1, x)
    print(x)
    print(y)
    plt.plot(x, y)  
    plt.show()

def graph_all(formula, x_range, k_1_range):
	print(pi1_d)
	x = np.array(x_range)
	for k_1 in k_1_range:
		print(k_1)
		y = formula(k_1,x)
		print(y)
		plt.plot(x,y)
	#plt.show()

def frange(start, stop, step):
    i = start
    while i < stop:
        yield i
        i += step

# pi1_init = float(input('What is pi1_init?   '))
# pi2_init = float(input('What is pi2_init?   '))
# xi1 = float(input('What is xi1?   '))
# xi2 = float(input('What is xi2?   '))
# xi1_d = float(input('What is xi1_delta?   '))
# xi2_d = float(input('What is xi2_delta?   '))
# r1 = float(input('What is region i?   '))
# r2 = float(input('What is region j?   '))
# r1_d = float(input('What is the attack target region l?   '))
# r2_d = float(input('What is the attack target region m?   '))

pi1 = .4
pi2 =  .4
pi1_d = 0.467
pi2_d = 0.467
xi1 =   .75
xi2 =   .75
xi1_d =   .75
xi2_d =   .75
r1 =   3
r2 =   8
r1_d =   3
r2_d =   8

print(pi1)
print(pi2)
print(pi1_d) 
print(pi2_d)
print(xi1)
print(xi2)
print(xi1_d)
print(xi2_d)
print(r1)
print(r2)
print(r1_d)
print(r2_d)

#default values are:
# pi1_init = float(14/30) = .467
# pi2_init = float(14/30) = .467
# lost time = 2 s 
# therefore, (30-2)/2 = 14s = effective green time
# therefore, effective green time ratio 1 and 2: pi1 = pi2 = 14/30
#xi1,2 init = 0.75
#jam density = 150vpmile
#k is 85
pi1_d = 0.1
pi2_d = 0.1
T = float(30) # seconds -> hours
t_update = float(0.01/3600) # updating time step; secs -> hours
#xi1 = .55
#xi2 = .55
k_j = 150 # jam density ranges from 185-250 per mile per lane typically; they used 60 in example
k_c = k_j / 5 # critical density is approx 1/5 kj
k_1_init = 50 # depends on initial region now
k_2_init = 100 # depends on initial region now
k_default = 85 # average overall density cars per mile for both rings (should be constant) should be k >= kj/2 to enable multivaluedness
# also depends on initial region
v_f = 70 # miles per hour ; free flow speed 
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
#pi1_ds = np.arange(pi1_init-1, 1-pi1_init , 0.01)
def dk1_dt_r1(k_1, t): 
	return (t * -1) * (1-xi1) * v_f * k_1 / L 
def dk1_dt_r2(k_1, t): 
	return (t * -1) * (1-xi1) * v_f * k_c / L 
def dk1_dt_r3 (k_1, t):
	return (t * -1) * (k_j-k_1) * ((1-xi1)*v_f*k_c) / (L*xi1*(k_j-k_1)) 
def dk1_dt_r4(k_1, t): 
	return (t * -1 ) * (v_f * k_c) / (L* (k_j-k_c)) * (k_j - 2*k+ k_1) 
def dk1_dt_r5(k_1, t): 
	return t * ((1-xi2) * v_f) / L * (2*k-k_1)  
def dk1_dt_r6 (k_1, t):
	return t * ((1-xi2) * v_f) / L * k_c
def dk1_dt_r7(k_1, t): 
	return t * ((1-xi2) * v_f * k_c) / (L* xi2 * (k_j-k_c)) * (k_j-2*k+k_1) 
def dk1_dt_r8(k_1, t): 
	return t * (v_f*k_c)/(L*(k_j-k_c)) * (k_j-k_1) 


g1 = (1-xi1)*v_f/L
g2 = ((1-xi1)*v_f*k_c)/(L*xi1*(k_j-k_c))
g3 = v_f*k_c / (L*(k_j-k_c))
g4 = (1-xi2)*v_f/L
g5 =  ((1-xi2)*v_f*k_c)/(L*xi2*(k_j-k_c))

# region intervals, mins and maxs
r_1_k_1 = interval()
r_2_k_1 = interval()
r_3_k_1 = interval()
r_4_k_1 = interval()
r_5_k_1 = interval()
r_6_k_1 = interval()
r_7_k_1 = interval()
r_8_k_1 = interval()
r_1_k = interval()
r_2_k = interval()
r_3_k = interval()
r_4_k = interval()
r_5_k = interval()
r_6_k = interval()
r_7_k = interval()
r_8_k = interval()

def def_regions(k_1, k):
	global r_1_k_1
	r_1_k_1 = interval([0,k_c])
	global r_1_k
	r_1_k = interval ([(k_1/2), (((1-xi1)*k_j - (2-xi1)*k_c) * k_1)/(2*k_c)])
	global r_2_k_1
	r_2_k_1 = interval([k_c,k_j-xi1*(k_j-k_c)])
	global r_2_k
	r_2_k = interval ([(k_1/2), (xi1*k_j + (1-xi1)*k_c + k_1)/2])
	global r_3_k_1	
	r_3_k_1 = interval([k_j-xi1*(k_j-k_c), k_j])
	global r_3_k
	r_3_k = interval ([(k_1/2), ((2*xi1-1)/(2*xi1) * k_j + (k_1/(2*xi1)))])
	global r_4_k_1
	r_4_k_1 = interval (0, k_j)
	global r_4_k
	r_4_k = interval([max((k_j/2 -(((1-xi1)*k_j - (2-xi1)*k_c)*k_1)/(2*k_c)), ((xi1*k_j) + ((1-xi1)*k_c) + k_1)/2, (2*xi1 - 1)/(2*xi1) + (k_1)/(2*xi1) ),(k_1+k_j)/2])
	global r_5_k_1
	r_5_k_1 = interval([0, k_j])
	global r_5_k
	r_5_k = interval ([(k_1/2), min(((k_1 + k_c)/2), ( k_1/2 + (k_c*(k_j-k_1))/(2*(1-xi2)*(k_j-k_c)) ))])
	global r_6_k_1
	r_6_k_1 = interval([0, k_j-(1-xi2)*(k_j-k_c)])
	global r_6_k
	r_6_k = interval ([((k_1+k_c)/2), ((k_j + k_1 - xi2*(k_j-k_c))/2) ])
	global r_7_k_1
	r_7_k_1 = interval([0,k_j])
	global r_7_k
	r_7_k = interval ([max(((k_j + k_1 - xi2*(k_j-k_c))/2), ((1-2*xi2)*k_j + k_1)/(2*(1-xi2))), (k_1+k_j)/2])
	global r_8_k_1
	r_8_k_1 = interval([k_j-(1-xi2)*(k_j-k_c), k_j])
	global r_8_k
	r_8_k = interval ([( k_1/2 + (k_c*(k_j-k_1))/(2*(1-xi2)*(k_j-k_c)) ) , ((1-2*xi2)*k_j + k_1)/(2*(1-xi2))])

# each of the # region intervals, mins and maxsfollowing must be altered according to pyinterval
# check if (k_1,k) is in the region by checking if k_1 satisfies R_i_k_1
# and check if (k_1,k) is in the region by checking if k satisfies R_i_k
# format: create intervals corresponding to min and max of region i
# since there are open intervals, this might not be doable; unless we do a precheck of the endpoints .
# return ( k_1 in r_k_1 and k in r_k) 
def region_1_check(k_1,k):
	in_region = 0
	if k_1 in r_1_k_1 and k in r_1_k and k_1 < r_1_k_1.extrema[1][0] and k_1 > r_1_k_1.extrema[0][0]:
		in_region = 1
	else:
		in_region = 0
	return (in_region)
#Aho

def region_2_check(k_1,k):
	in_region = 0
	if k_1 in r_2_k_1 and k in r_2_k  and  k_1 < r_2_k_1.extrema[1][0]:
		in_region = 1
	else:
		in_region = 0
	return (in_region)
#Aho
def region_3_check (k_1,k):
	in_region = 0
	if k_1 in r_3_k_1 and k in r_3_k and k_1 < r_3_k_1.extrema[1][0] and k_1 > r_3_k_1.extrema[0][0]:
		in_region = 1
	else:
		in_region = 0
	return (in_region)
def region_4_check (k_1,k):
	in_region = 0
	if k_1 in r_4_k_1 and k in r_4_k and k_1 < r_4_k_1.extrema[1][0] and k_1 > r_4_k_1.extrema[0][0] and k > r_4_k.extrema[0][0]:
		in_region = 1
	else:
		in_region = 0
	return (in_region)
def region_5_check (k_1,k):
	in_region = 0
	if k_1 in r_5_k_1 and k in r_5_k:
		in_region = 1
	else:
		in_region = 0
	return (in_region)
def region_6_check (k_1,k):
	in_region = 0
	if k_1 in r_6_k_1 and k in r_6_k and k > r_6_k.extrema[0][0]:
		in_region = 1
	else:
		in_region = 0
	return (in_region)
def region_7_check(k_1,k):
	in_region = 0
	if k_1 in r_7_k_1 and k in r_7_k and k > r_7_k.extrema[0][0]:
		in_region = 1
	else:
		in_region = 0
	return (in_region)
def region_8_check (k_1,k):
	in_region = 0
	if k_1 in r_8_k_1 and k in r_8_k and k_1 > r_8_k_1.extrema[0][0] and k > r_8_k.extrema[0][0] and k < r_8_k.extrema[1][0]:
		in_region = 1
	else:
		in_region = 0
	return (in_region)

# Green Time Ratio Modification Attack:

#given that initial state is 
#1) same green time delta pi1 delta = pi2 delta -> same effect as changing cycle time for a given fixed state
poinc_map_3_8 = k_j * (1-np.exp((g2-g3)*pi1*T/(3600))) + k_1_init*np.exp((g2-g3)*pi1*T/(3600))

#graph poinc map 3_8 with normal pi
# append new k_1 to end of the 
print(g2-g3)
def poinc_map_3_8_delta (k_1, iters):
	return (k_j * (1-np.exp((g2-g3)*pi1_d*T/(3600)*iters)) + k_1*np.exp((g2-g3)*pi1_d*T/(3600)*iters))

#y = poinc_map_3_8_delta(130, np.array(range(0,40)))
#print(y)
#graph(poinc_map_3_8_delta, range(0,40), 130)
step = 2.0
# # initial density should be the same as the minimum value for 
def_regions(k_1_init, k_default)
print(k_1_init)
print(k_default)
# # we need to start in the (2,6) region
# # we can use k_1 min, k_1 max if possible 
# # and accordingly adopt all valid k values within this region (R_2_k intersect R_6_k)


k_1_min = r_3_k_1.extrema[0][0]
k_1_max = r_3_k_1.extrema[1][0]
# for the above sets, we need to determine which values of k_1 satisfy the DEO (2,6)
# we can use the dk1/dt formulas and region checks to determine this


# print(k_1_min)
# print(k_1_max)
k_1_vals = frange(130,k_1_max, step)
k_1_vals = [130]
pi_vals = frange(0.2, 5.0, 0.4)
cycle_length_vals = frange(30, 60, 5)
cycle_range = range(0,40)
# for val in pi_vals:
# 	pi1_d = val
# 	graph_all(poinc_map_3_8_delta, cycle_range, k_1_vals)
# 	#plt.show()
# plt.show()

pi1_d = 0.467
for val in cycle_length_vals:
	T= val
	graph_all(poinc_map_3_8_delta, cycle_range, k_1_vals)
	#plt.show()
plt.show()
