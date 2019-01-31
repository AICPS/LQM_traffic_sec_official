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


def graph(formula, x_range, k_1):  
    x = np.array(x_range)  
    y = formula(k_1, x)
    print(x)
    print(y)
    plt.plot(x, y)
    plt.show()

def graph_all(formula, x_range, attack_start, attack_end, k_1_range, all_data_ws, col ):
	# name = 'pi' + str(pi1_d)
	# ws = wb.add_worksheet(name)
	# print(pi1_d)
	x = np.array(x_range)
	i = 0
	for k_1 in k_1_range:
		print(k_1)
		y = formula(k_1,x)
		# print(y)
		plt.plot(x,y)

		for row, data in enumerate(y):
			all_data_ws.write(row, col, data)
		
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

# settings here for initial settings/ state and attack modifications
pi1 = .35
pi2 =  .35
xi1 =   .75
xi2 =   .75
pi1_d = 0.467
pi2_d = 0.467
xi1_d =   .75
xi2_d =   .75
pi = pi1

print(pi1)
print(pi2)
print(pi1_d) 
print(pi2_d)
print(xi1)
print(xi2)
print(xi1_d)
print(xi2_d)

#default values are:
# pi1_init = float(14/30) = .467
# pi2_init = float(14/30) = .467
# lost time = 2 s 
# therefore, (30-2)/2 = 14s = effective green time
# therefore, effective green time ratio 1 and 2: pi1 = pi2 = 14/30
#xi1,2 init = 0.75
#jam density = 150vpmile
#k is 85
# constraints
#pi1 + pi2 + offset = 1
#xi1, xi2 > 0.5

T = float(30) # seconds -> hours
t_update = float(0.01) # updating time step; secs -> hours
#xi1 = .55
#xi2 = .55
k_j = 150 # jam density ranges from 185-250 per mile per lane typically; they used 60 in example
k_c = k_j / 5 # critical density is approx 1/5 kj
k_1_init = 130 # depends on initial region now
k_default = 85 # average overall density cars per mile for both rings (should be constant) should be k >= kj/2 to enable multivaluedness
k_2_init = 2*k_default - k_1_init # depends on initial region now
k = k_default
# also depends on initial region
v_f = 60 # miles per hour ; free flow speed 
L = 1 # mile length of each ring in miles
total_cycles = 70

#attack model information
attack_cycle_start = 5
attack_cycle_end = total_cycles
attack_reg = (3,8)

g1 = (1-xi1)*v_f/L
g2 = ((1-xi1)*v_f*k_c)/(L*xi1*(k_j-k_c))
g3 = v_f*k_c / (L*(k_j-k_c))
g4 = (1-xi2)*v_f/L
g5 =  ((1-xi2)*v_f*k_c)/(L*xi2*(k_j-k_c))


def error_function ():
	#error probability distribution function
	#sample value
	return .1

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
	if k_1 in r_3_k_1 and k in r_3_k and k_1 <= r_3_k_1.extrema[1][0] and k_1 => r_3_k_1.extrema[0][0]:
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

# given current values of k1 and k, return the current orbit based on the region boundaries
# and based on dk1/dt
def return_deo(k_1,k):
	deo = (-1,-1)
	i = -1
	j = -1
	i = curr_region(k_1,k,1)
	#print(update_k1(k_1,k,i)) 
	k_1 = max(0, update_k1(k_1,k,i)) # update for phase 1; region i
	j = curr_region(k_1,k,2)
	k_1 = min(k_j, update_k1(k_1,k,j)) # udpate for phase 2; region j
	deo = (i,j)
	return deo

k1_max  = k_j - (xi1)*(k_j-k_c) # for it to be in region 2
k1_min = k_c
k_1 = k_j - (xi1)*(k_j-k_c) - 1
k = (xi1 * k_j + (2-xi1)*k_c) / 2
def_regions(k_1, k)
k_values = r_2_k & r_6_k # get the values for k that satisfy the intersection of these intervals

#given that initial state is 
#1) same green time delta pi1 delta = pi2 delta -> same effect as changing cycle time for a given fixed state

# focus on asymptotic stationary states first then lyapunov, then unstable
# asym when xi > 0.5; lyap when xi = 0.5; unstable when xi < 0.5
poinc_map_3_8 = k_j * (1-np.exp((g2-g3)*pi*T/(3600))) + k_1_init*np.exp((g2-g3)*pi*T/(3600))
# asym when xi > 0.5; lyap when xi = 0.5; unstable when xi < 0.5
poinc_map_4_7 = (k_j - 2*k) *(np.exp(g5*pi*T/3600)-1) + k_1_init*np.exp((g5-g1)*pi*T/(3600))
# asym when xi > 0.5; 
poinc_map_3_5 = k_j * (1-np.exp(g2*pi*T/3600)) * np.exp(-g4*pi*T/3600) + 2*k*(1-np.exp(-g4*pi*T/3600)) + (k_1_init * np.exp((g2-g4)*pi*T/3600))
# asym when xi > 0.5; 
poinc_map_1_7 = (k_j - 2*k)*(np.exp(g5*pi*T/3600)-1) + k_1_init * np.exp((g5-g3)*pi*T/3600)
# asym when xi < 0.5; 
poinc_map_4_8 = k_j * (np.exp(-2*g3*pi*T/3600)-2*np.exp(-g3*pi*T/3600) + 1) - 2*k *(np.exp(-2*g3*pi*T/3600) -np.exp(-g3*pi*T/3600)) + k_1_init * np.exp(-2*g3*pi*T/3600)
# unstable when xi > 0.5
poinc_map_3_7 = k_j * (2*np.exp(g2*pi*T/3600) - np.exp(2*g2*pi*T/3600) -1) - 2*k*(np.exp(g2*pi*T/3600) - 1) + k_1_init *np.exp(2*g2*pi*T/3600)

#graph poinc map 3_8 with normal pi
# append new k_1 to end of the 

# asym when xi > 0.5; lyap when xi = 0.5; unstable when xi < 0.5
def poinc_map_3_8_delta (k_1, n, pi_arg):
	return (k_j * (1-np.exp((g2-g3)*pi_arg*T/(3600)*n)) + k_1*np.exp((g2-g3)*pi_arg*T/(3600)*n))
# asym when xi > 0.5; lyap when xi = 0.5; unstable when xi < 0.5
def poinc_map_4_7_delta (k_1, n, pi_arg):
	return (k_j - 2*k) *(np.exp(g5*pi_arg*T/3600*n)-1) + k_1*np.exp((g5-g1)*pi_arg*T/(3600)*n)
# asym when xi > 0.5; 
def poinc_map_3_5_delta (k_1, n, pi_arg):
	return k_j * (1-np.exp(g2*pi_arg*T/3600*n)) * np.exp(-g4*pi_arg*T/3600*n) + 2*k*(1-np.exp(-g4*pi_arg*T/3600*n)) + (k_1 * np.exp((g2-g4)*pi_arg*T/3600*n))
# asym when xi > 0.5; 
def poinc_map_1_7_delta (k_1, n, pi_arg):
	return (k_j - 2*k)*(np.exp(g5*pi_arg*T/3600*n)-1) + k_1 * np.exp((g5-g3)*pi_arg*T/3600*n)
# asym when xi < 0.5; 
def poinc_map_4_8_delta (k_1, n, pi_arg):
	return k_j * (np.exp(-2*g3*pi_arg*T/3600*n)-2*np.exp(-g3*pi_arg*T/3600*n) + 1) - 2*k *(np.exp(-2*g3*pi_arg*T/3600*n) -np.exp(-g3*pi_arg*T/3600*n)) + k_1 * np.exp(-2*g3*pi_arg*T/3600*n)
# unstable when xi > 0.5
def poinc_map_3_7_delta (k_1, n, pi_arg):
	return k_j * (2*np.exp(g2*pi_arg*T/3600*n) - np.exp(2*g2*pi_arg*T/3600*n) -1) - 2*k*(np.exp(g2*pi_arg*T/3600*n) - 1) + k_1 *np.exp(2*g2*pi_arg*T/3600*n)

def compute_k1_array(k_1, R, n_range, pi_arg):
	print('r is: ', R)
	k1_arr = []
	if R == (3,8):
		for n in n_range:
			k1 = poinc_map_3_8_delta (k_1, n, pi_arg)
			if(k1 < 0):
				k1 = 0
			elif(k1 > k_j):
				k1 = k_j
			k1_arr.append(k1)
		return k1_arr
	elif R==(4,7): 
		for n in n_range:
			k1 = poinc_map_4_7_delta (k_1, n, pi_arg)
			if(k1 < 0):
				k1 = 0
			elif(k1 > k_j):
				k1 = k_j
			k1_arr.append(k1)
		return k1_arr
	elif R == (3,5):
		for n in n_range:
			k1 = poinc_map_3_5_delta (k_1, n, pi_arg)
			if(k1 < 0):
				k1 = 0
			elif(k1 > k_j):
				k1 = k_j
			k1_arr.append(k1)
		return k1_arr
	elif R == (1,7):
		for n in n_range:
			k1 = poinc_map_1_7_delta (k_1, n, pi_arg)
			if(k1 < 0):
				k1 = 0
			elif(k1 > k_j):
				k1 = k_j
			k1_arr.append(k1)
		return k1_arr
	elif R == (4,8):
		for n in n_range:
			k1 = poinc_map_4_8_delta (k_1, n, pi_arg)
			if(k1 < 0):
				k1 = 0
			elif(k1 > k_j):
				k1 = k_j
			k1_arr.append(k1)
		return k1_arr
	elif R == (3,7):
		for n in n_range:
			k1 = poinc_map_3_7_delta (k_1, n, pi_arg)
			if(k1 < 0):
				k1 = 0
			elif(k1 > k_j):
				k1 = k_j
			k1_arr.append(k1)
		return k1_arr
	else:
		print('Invalid DEO')
		return -1.0

#y = poinc_map_3_8_delta(130, np.array(range(0,40)))
#print(y)
#graph(poinc_map_3_8_delta, range(0,40), 130)
#step = 1.0
# # initial density should be the same as the minimum value for 
def_regions(k_1_init, k_default)
k_1_min = r_3_k_1.extrema[0][0]
k_1_max = r_3_k_1.extrema[1][0]
# print(k_1_min)
# print(k_1_max)
#k_1_vals = frange(130,k_1_max, step)
k_1_vals = [130]
k_1_init = 130
k = 85
#pi_vals = frange(0.2, 0.4, 0.4)
attack_region = (3,8)

timestep = 1
time_range = range(0, total_cycles, 1) # number of simulation cycles
pi_vals = frange(0.44, 0.52, 0.02)
#t_vals = frange(20, 60, 5)
# iterate through each modified cycle length
#for val in t_vals:
file_name = 'D:/Anthony/Dropbox/Anthony_Research/Automotive Security/Writing/NDSS-2019/Data/same_GTR_3_8_poinc_map_initpi_'+str(pi)+'_xi_'+str(xi1)+'.xlsx'
wb = xlsw.Workbook(file_name)
#wb = xlsw.Workbook('./data/alldata_same_gtr.xlsx')
all_data_ws = wb.add_worksheet('data')
pi_val_sheet = wb.add_worksheet('pi_vals')
for row, data in enumerate(time_range):
	print(data)
	all_data_ws.write(row, 0, data)
col = 1

# compute values for initial pi value
pi_val_sheet.write(0, 0, pi)
x = list(time_range)
y = compute_k1_array(k_1_init, attack_region, range(0,total_cycles, timestep), pi)
#print('y: ', y)
print('length of time range ' , len(x))
# print('len(x) ', len(x))
# print('len(y) ', len(y))

plt.plot(x,y)
for row, data in enumerate(y):
	all_data_ws.write(row, col, data)
col += 1
# for each value of pi that we are testing
for val in pi_vals:
	y1 = []
	y_attack = []
	y2 = []
	y = []
	#print('T_d', val)
	print('pi_d: ', val)
	pi_val_sheet.write(0, col, val)
	# pi1_d = val 
	# pi2_d = val
	#T = val
	# graph_all(poinc_map_3_8_delta, cycle_range, k_1_vals, all_data_ws, col)
	# name = 'pi' + str(pi1_d)
	# ws = wb.add_worksheet(name)
	# print(pi1_d)

	#TODO: CHANGE THIS PART TO INCLUDE ATTACK AT START OF ATTACK AND END OF ATTACK

	y1 = compute_k1_array(k_1_init, attack_region,range(0,attack_cycle_start, timestep), pi)
	y_attack = compute_k1_array(k_1_init, attack_region, range(attack_cycle_start,attack_cycle_end, timestep), val)
	y2 = compute_k1_array(k_1_init, attack_region, range(attack_cycle_end,total_cycles, timestep), pi)
	y = y1 + y_attack + y2
	#print('y: ', y)
	plt.plot(x,y)
	for row, data in enumerate(y):
		all_data_ws.write(row, col, data)

	# for k_1 in k_1_range:
	# 	print(k_1)
	# 	y = formula(k_1,x)
	# 	# print(y)
	# 	plt.plot(x,y)

	# 	for row, data in enumerate(y):
	# 		all_data_ws.write(row, col, data)
	#plt.show()
	col +=1 

plt.show()
wb.close()
