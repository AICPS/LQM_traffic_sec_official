import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
from scipy.integrate import odeint
import sys 
from interval import interval, inf, imath, fpu
from scipy.integrate import odeint
import xlwt
from datetime import datetime
import xlsxwriter as xlsw

# def f (y, t, params):
# 	k1, k2 = y
# 	k, L, xi1, xi2, delta1, delta2, g1, g2, n, T = params
# 	#print("k1 ",k1, "at t ", t)
# 	#print("k2 ", k2, " at t ", t)
# 	derivs = [-(1-xi1*g1(k1,k2,t,n,T))/L + (1-xi2*g2(k1,k2,t,n,T))/L, -(1-xi2*g2(k1,k2,t,n,T))/L + (1-xi1*g1(k1,k2,t,n,T))/L] # out flow - in flow
# 	return derivs

def frange(start, stop, step):
    i = start
    while i < stop:
        yield i
        i += step

# settings here for initial settings/ state and attack modifications
pi1 = .35
pi2 =  .35
xi1 =   .75
xi2 =   .75
pi = pi1
pi1_d = 0.6
pi2_d = 0.4
pi_d_arr = [pi]
for elem in list(frange(0.44, 0.52, 0.02)):
	pi_d_arr.append(elem)
#pi_d_arr = frange(0.30, 0.48, 0.2)
#pi2_d_arr = frange(0.30, 0.48, 0.2)
xi1_d =   .75
xi2_d =   .75

#attack model information
attack_cycle_start = 5
attack_cycle_end = 99
attack_reg = (3,8)

tInc = 0.01
phase_offset = 0.0 

T = float(60) # seconds -> hours
t_update = float(0.01) # updating time step; secs -> hours
#xi1 = .55
#xi2 = .55
k_j = 150 # jam density ranges from 185-250 per mile per lane typically; they used 60 in example
k_c = k_j / 5 # critical density is approx 1/5 kj
k_default = 85 # average overall density cars per mile for both rings (should be constant) should be k >= kj/2 to enable multivaluedness
k = k_default
v_f = 60 # miles per hour ; free flow speed 
L = 1 # mile length of each ring in miles
k1_min = 130
k1_max = 150
w_a = 40 #miles per hour

#gammas
g1 = (1-xi1)*v_f/L
g2 = ((1-xi1)*v_f*k_c)/(L*xi1*(k_j-k_c))
g3 = v_f*k_c / (L*(k_j-k_c))
g4 = (1-xi2)*v_f/L
g5 =  ((1-xi2)*v_f*k_c)/(L*xi2*(k_j-k_c))

def dk1_dt_r1(k_1,t): 
	#return -1 * (1-xi1) * v_f * k_1 / L * pi1_d * T * t/(3600)
	return (-1 * g1 * k_1) * pi1_d * t/(3600)
def dk1_dt_r2(k_1,t): 
	#return -1 * (1-xi1) * v_f * k_c / L * pi1_d * T * t/(3600)
	return (-1 * g1 * k_c) * pi1_d * t/(3600)
def dk1_dt_r3(k_1,t):
	#return -1 * (k_j-k_1) * ((1-xi1)*v_f*k_c) / (L*xi1*(k_j-k_1)) * pi1_d * T * t/(3600)
	return (g2*k_1 - g2*k_j) * pi1_d * t/(3600)
def dk1_dt_r4(k_1,t): 
	#return -1 * (v_f * k_c) / (L* (k_j-k_c)) * (k_j - 2*k+ k_1) * pi1_d * T * t/(3600)
	return ((-1 * g3 * k_1) - (g3 * (k_j - 2 * k))) * pi1_d * t/(3600)
def dk1_dt_r5(k_1,t): 
	#return ((1-xi2) * v_f) / L * (2*k-k_1)  * pi2_d * T * t/(3600)
	return ((- 1 * g4 * k_1) + (2 * g4 * k)) * pi2_d * t/(3600)
def dk1_dt_r6(k_1,t):
	#return ((1-xi2) * v_f * k_c) / L * pi2_d * T * t/(3600)
	return (g4 * k_c) * pi2_d * t/(3600)
def dk1_dt_r7(k_1,t): 
	#return  ((1-xi2) * v_f * k_c) / (L* xi2 * (k_j-k_c)) * (k_j-2*k+k_1) * pi2_d * T * t/(3600)
	return ((g5 * k_1) + g5 * (k_j - 2 * k)) * pi2_d  * t/(3600)
def dk1_dt_r8(k_1,t): 
	#return  (v_f*k_c)/(L*(k_j-k_c)) * (k_j-k_1) * pi2_d * T * t/(3600)
	return ((-1 * g3 * k_1) + (g3 * k_j)) * pi2_d * t/(3600)


def dk1_dt_r2_6(k_1,t): 
	return dk1_dt_r2(k_1,t) + dk1_dt_r6(k_1,t)
def dk1_dt_r4_7(k_1,t): 
	return dk1_dt_r4(k_1,t) + dk1_dt_r7(k_1,t)
def dk1_dt_r3_6(k_1,t):
	return dk1_dt_r3(k_1,t) + dk1_dt_r6(k_1,t)
def dk1_dt_r3_7(k_1,t): 
	return dk1_dt_r3(k_1,t) + dk1_dt_r7(k_1,t)
def dk1_dt_r3_8(k_1,t): 
	return dk1_dt_r3(k_1,t) + dk1_dt_r8(k_1,t)
def dk1_dt_r4_8(k_1,t):
	return dk1_dt_r4(k_1,t) + dk1_dt_r8(k_1,t)


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
	r_4_k_1 = interval ([0, k_j])
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
	r_8_k = interval ([(k_1/2) - (k_c*(k_j-k_1))/(2*(1-xi2)*(k_j-k_c)), ((1-2*xi2)*k_j + k_1)/(2*(1-xi2))])


def update_k1(k_1,k,R1, t):
	print('R1 is ', R1)
	def_regions(k_1, k)
	k0 = k_1
	if R1 == 1:
		k1_r1 = odeint(dk1_dt_r1,k0,t)
		k_1 = k1_r1[-1]
		return k_1
	elif R1 == 2:
		k1_r2 = odeint(dk1_dt_r2,k0,t)
		k_1 = k1_r2[-1]
		return k_1
	elif R1 == 3:
		print('In r1 block')
		k1_r3 = odeint(dk1_dt_r3,k0,t)
		print('k1_r3 ', k1_r3)
		k_1 = k1_r3[-1]
		return k_1
	elif R1 == 4:
		k1_r4 = odeint(dk1_dt_r4,k0,t)
		k_1 = k1_r4[-1]
		return k_1
	elif R1 == 5:
		k1_r5 = odeint(dk1_dt_r5,k0,t)
		k_1 = k1_r5[-1]
		return k_1
	elif R1 == 6:
		k1_r6 = odeint(dk1_dt_r6,k0,t)
		k_1 = k1_r6[-1]
		return k_1
	elif R1 == 7:
		k1_r7 = odeint(dk1_dt_r7,k0,t)
		k_1 = k1_r7[-1]
		return k_1
	elif R1 == 8:
		k1_r8 = odeint(dk1_dt_r8,k0,t)
		k_1 = k1_r8[-1]
		return k_1	
	else:
		print('error in update_k1')
		return -1.0

def update_k1_deo(k_1,k,R1,R2):
	def_regions(k_1, k)
	t = np.array(range(0,2))
	#k0 = k_1
	k1 = update_k1(k_1,k,R1)
	k1 = update_k1(k_1,k,R2)
	return k1

#TODO	
def Q(k_a): # flow-density relationship
	return min(v_f*k_a, w_a*(k_j-k_a))
	
def D(k_a): # demand - aka out-ward flow (into junction)
	demand = 0
	if(k_a >= 0 and k_a <= k_c):
		demand = Q(k_a)
	elif (k_a > k_c and k_a <= k_j): # k_a is between k_c and k_j
		demand = Q(k_c)
	else:
		#print("")
		demand = 0
	return demand

def S(k_a): # supply - aka in-ward flow (into link from junction)
	supply = 0
	if(k_a >= 0 and k_a <= k_c):
		supply = Q(k_c)
	elif (k_a > k_c and k_a <= k_j): # k_a is between k_c and k_j
		supply = Q(k_a)
	else:
		#print("")
		demand = 0
	return supply

# do i create an array for all values of g1

def outflow1(k1,k2):
	val = min(D(k1), S(k1)/xi1, S(k2)/(1-xi1))
	#print("g1 at time", t, " is ", val)
	return val

def outflow2(k1,k2):
	val = min(D(k2), S(k2)/xi2, S(k1)/(1-xi2))
	#print("g2 at time", t, " is ", val)
	return val

#def integrand (t, flow1, flow2):

# def compute_asym_flow(lower, upper, k1_1, k2_1, k1_2, k2_2, T):
# 	answer1 = quad(outflow1, lower, upper, args=(k1_1,k2_1,T))
# 	answer2 = quad(outflow2, lower, upper, args=(k1_2,k2_2,T))
# 	flux1 = answer1[0]
# 	flux2 = answer2[0]
# 	q = (flux1+flux2)/(2*T) # the average network flow during cycle
# 	return q
	#second element is the error
	#integral of g1 from t-T to t + integral of g2 from t-T to t / 2*T 

def compute_flow(p, k1,k2):
	if(p == 1):
		return outflow1(k1,k2)
	else:
		return outflow2(k1,k2)

def error_function ():
	#error probability distribution function
	#sample value
	return .1

# each of the # region intervals, mins and maxsfollowing must be altered according to pyinterval
# check if (k_1,k) is in the region by checking if k_1 satisfies R_i_k_1
# and check if (k_1,k) is in the region by checking if k satisfies R_i_k
# format: create intervals corresponding to min and max of region i
# since there are open intervals, this might not be doable; unless we do a precheck of the endpoints .
# return ( k_1 in r_k_1 and k in r_k) 
def region_1_check(k_1,k):
	def_regions(k_1, k)
	in_region = 0
	if k_1 in r_1_k_1 and k in r_1_k and k_1 < r_1_k_1.extrema[1][0] and k_1 > r_1_k_1.extrema[0][0]:
		in_region = 1
	else:
		in_region = 0
	return (in_region)
#Aho

def region_2_check(k_1,k):
	def_regions(k_1, k)
	in_region = 0
	if k_1 in r_2_k_1 and k in r_2_k  and  k_1 < r_2_k_1.extrema[1][0]:
		in_region = 1
	else:
		in_region = 0
	return (in_region)
#Aho
def region_3_check (k_1,k):
	def_regions(k_1, k)
	in_region = 0
	if k_1 in r_3_k_1 and k in r_3_k and k_1 < r_3_k_1.extrema[1][0] and k_1 > r_3_k_1.extrema[0][0]:
		in_region = 1
	else:
		in_region = 0
	return (in_region)

def region_4_check (k_1,k):
	def_regions(k_1, k)
	in_region = 0
	if k_1 in r_4_k_1 and k in r_4_k and k_1 < r_4_k_1.extrema[1][0] and k_1 > r_4_k_1.extrema[0][0] and k > r_4_k.extrema[0][0]:
		in_region = 1
	else:
		in_region = 0
	return (in_region)

def region_5_check (k_1,k):
	def_regions(k_1, k)
	in_region = 0
	if k_1 in r_5_k_1 and k in r_5_k:
		in_region = 1
	else:
		in_region = 0
	return (in_region)

def region_6_check (k_1,k):
	def_regions(k_1, k)
	in_region = 0
	if k_1 in r_6_k_1 and k in r_6_k and k > r_6_k.extrema[0][0]:
		in_region = 1
	else:
		in_region = 0
	return (in_region)

def region_7_check(k_1,k):
	def_regions(k_1, k)
	in_region = 0
	if k_1 in r_7_k_1 and k in r_7_k and k > r_7_k.extrema[0][0]:
		in_region = 1
	else:
		in_region = 0
	return (in_region)

def region_8_check (k_1,k):
	def_regions(k_1, k)
	in_region = 0
	if k_1 in r_8_k_1 and k in r_8_k and k_1 > r_8_k_1.extrema[0][0] and k > r_8_k.extrema[0][0] and k < r_8_k.extrema[1][0]:
		in_region = 1
	else:
		in_region = 0
	return (in_region)

def curr_region(k_1,k,p):
	def_regions(k_1, k)
	if p == 1:
		if(region_1_check(k_1,k)):
			return 1
		elif(region_2_check(k_1,k)):
			return 2
		elif(region_3_check(k_1,k)):
			return 3
		elif(region_4_check(k_1,k)):
			return 4
		else:
			return -1
	elif p == 2:
		if(region_5_check(k_1,k)):
			return 5
		elif(region_6_check(k_1,k)):
			return 6
		elif(region_7_check(k_1,k)):
			return 7
		elif(region_8_check(k_1,k)):
			return 8
		else:
			return -1
	else:
		print('error in curr_region')
		return -1

# plotting setup
#if(len(cycle_range) == len(k1_arr)): 
fig = plt.figure(1,figsize=(8, 8))
# density vs time
ax1 = fig.add_subplot(311)
# flow vs time
ax2 = fig.add_subplot(312)
# flow vs density - fundamental diagram
ax3 = fig.add_subplot(313)

#XLS writing setup
file_name = 'D:/Anthony/Dropbox/Anthony_Research/Automotive Security/Writing/NDSS-2019/Data/same_GTR_3_8_initpi_'+str(pi)+'_xi_'+str(xi1)+'.xlsx'
wb = xlsw.Workbook(file_name)
k1_data_ws = wb.add_worksheet('data_k1')
k2_data_ws = wb.add_worksheet('data_k2')
q_data_ws = wb.add_worksheet('data_q')
pi_val_sheet = wb.add_worksheet('pi_vals')

# range of cycles
cycle_range = range(1,30)
for row, data in enumerate(cycle_range):
	k1_data_ws.write(row, 0, data)
	k2_data_ws.write(row, 0, data)
	q_data_ws.write(row, 0, data)

# for col, data in enumerate(pi_d_arr):
# 	pi_val_sheet.write(0, col, data)

q_arr = [] # array of flows during each time step
k1_arr = [] # array of average link density of k1 during each time step
k2_arr = [] # array of average link density of k2 during each time step
in_attack_reg = 0
# for different init densities
for k1_0 in range(k1_min,k1_max+1,2):
	q_arr = [] # array of flows during each time step
	k1_arr = [] # array of average link density of k1 during each time step
	k2_arr = []
	k2_0 = 2*k - k1_0
	k1 = k1_0
	k2 = k2_0
	print("k1_0 is: ", k1_0)
	print("k2_0 is: ", k2_0)
	print("k is: ", k)
	print("T is: ", T)
	#for each cycle, simulate each phase
	print('Length of pi_d_arr is: ', len(pi_d_arr))
	print('first region is: ', curr_region(k1,k,1))
	print('second region is: ', curr_region(k1,k,2))
	if(curr_region(k1,k,1) != attack_reg[0] or curr_region(k1,k,2) != attack_reg[1]):
		in_attack_reg = 0
	else:
		print('We found initial state in attack region')
		in_attack_reg = 1

	if(in_attack_reg == 1):
		col = 1
		for pi_d in pi_d_arr:
			q = compute_flow(1, k1_0, k2_0)
			q_arr.append(q)
			k1_arr.append(k1_0)
			k2_arr.append(k2_0)
			for cycle in cycle_range:
				print("Current cycle number is: ", cycle)
				for p in range(1,3):
					# Make time array
					print('p: ', p)
					curr_reg = curr_region(k1,k,p)
					print('Current region for these densities is: ', curr_reg)
					# no attack
					if(cycle <= attack_cycle_start):
						if(p == 1):
							t = list(frange(0,pi1*(T-phase_offset), tInc))
						else:
							t = list(frange(0,pi2*(T-phase_offset), tInc))
					elif(cycle > attack_cycle_end):
						if (cycle == attack_cycle_start):
							print('Stopping attack phase')
						if(p == 1):
							t = list(frange(0,pi1*(T-phase_offset), tInc))
						else:
							t = list(frange(0,pi2*(T-phase_offset), tInc))
					else:
						# attack (modified green time ratios)
						if (cycle == attack_cycle_start):
							print('We are in attack phase')
						if(p == 1):
							t = list(frange(0,pi_d*(T-phase_offset), tInc))
						else:
							t = list(frange(0,pi_d*(T-phase_offset), tInc))
					# for i in t:
					# 	print(i)
					# Call the ODE Solver
					print(k1)
					print(k2)
					#print(q)
					q = compute_flow(p, k1, k2) 
					k1 = update_k1(k1, k, curr_reg, t)
					if(k1 < 0):
						k1 = 0
					if(k1 > k_j):
						k1 = k_j
					k2 = 2*k - k1
					print(k1)
					print(k2)
					print(q)
					if(p == 1):
						q_arr.append(q)
					if(p == 1):
						k1_arr.append(k1)
					else:
						k2_arr.append(k2)
				# if(in_attack_reg == 0):
				#  	break
			# if(in_attack_reg == 1):
			#cycle_range[0:len(k1_arr)]	
			#cycle_range[0:len(k2_arr)]
			#cycle_range[0:len(q_arr)]		
			ax1.plot(range(0, len(k1_arr)), k1_arr, label="%.3f" % pi_d)
			ax2.plot(range(0, len(k2_arr)), k2_arr, label="%.3f" % pi_d)
			ax3.plot(range(0, len(q_arr)), q_arr, label="%.3f" % pi_d)
			
			for row, val in enumerate(k1_arr):
				k1_data_ws.write(row, col, val)
			for row, val in enumerate(k2_arr):
				k2_data_ws.write(row, col, val)
			for row, val in enumerate(q_arr):
				q_data_ws.write(row, col, val)
			col +=1 
			q_arr = [] # array of flows during each time step
			k1_arr = [] # array of average link density of k1 during each time step
			k2_arr = []
			k1 = k1_0
			k2 = k2_0

	# if(in_attack_reg == 1):
		# if(in_attack_reg == 1):
	print('We just finished simulation in attack region', attack_reg)
		# break
		# break;
	if(in_attack_reg == 1):
		break

		# else:
		# 	compute_asym_flow(lower, upper, k1_1, k2_1, k1_2, k2_2, T):

#Plot the results (if valid)
print('cycle range length: ', len(cycle_range))
print('k1 array length: ', len(k1_arr))
print('Length of pi_d_arr is: ', len(pi_d_arr))
ax1.set_xlabel('time t')
ax1.set_ylabel('density k1')
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels)
ax2.set_xlabel('time t')
ax2.set_ylabel('density k2')
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles, labels)
ax3.set_xlabel('time t')
ax3.set_ylabel('flow q')
handles, labels = ax3.get_legend_handles_labels()
ax3.legend(handles, labels)
plt.tight_layout()
plt.show()
wb.close()
# #if(len(cycle_range) == len(k1_arr)): 
# fig = plt.figure(1,figsize=(8, 8))
# # density vs time
# ax1 = fig.add_subplot(311)
# ax1.plot(cycle_range[0:len(k1_arr)], k1_arr)
# ax1.set_xlabel('time t')
# ax1.set_ylabel('density k1')
# # flow vs time
# ax2 = fig.add_subplot(312)
# ax2.plot(cycle_range[0:len(k2_arr)], k2_arr)
# ax2.set_xlabel('time t')
# ax2.set_ylabel('density k2')
# # flow vs density - fundamental diagram
# ax3 = fig.add_subplot(313)
# ax3.plot(cycle_range[0:len(q_arr)], q_arr)
# ax3.set_xlabel('time 2')
# ax3.set_ylabel('flow q')
# plt.tight_layout()
# plt.show()