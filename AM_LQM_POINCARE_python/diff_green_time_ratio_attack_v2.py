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


def beautiful_plot(range,data ):
	plt.plot(range, data)  
	# You typically want your plot to be ~1.33x wider than tall. This plot is a rare    
	# exception because of the number of lines being plotted on it.    
	# Common sizes: (10, 7.5) and (12, 9)    
	plt.figure(figsize=(12, 14))    
	  
	# Remove the plot frame lines. They are unnecessary chartjunk.    
	ax = plt.subplot(111)    
	ax.spines["top"].set_visible(False)    
	ax.spines["bottom"].set_visible(False)    
	ax.spines["right"].set_visible(False)    
	ax.spines["left"].set_visible(False)    
	  
	# Ensure that the axis ticks only show up on the bottom and left of the plot.    
	# Ticks on the right and top of the plot are generally unnecessary chartjunk.    
	ax.get_xaxis().tick_bottom()    
	ax.get_yaxis().tick_left()    

	plt.savefig("diff_green_time_attack.png", bbox_inches="tight") 
	plt.show() 


# pi1_init = float(input('What is pi1_init?   '))
# pi2_init = float(input('What is pi2_init?   '))
# xi1_d = float(input('What is xi1_d?   '))
# xi2_d = float(input('What is xi2_d?   '))
# xi1_d_d = float(input('What is xi1_d_delta?   '))
# xi2_d_d = float(input('What is xi2_d_delta?   '))
# r1 = float(input('What is region i?   '))
# r2 = float(input('What is region j?   '))
# r1_d = float(input('What is the attack target region l?   '))
# r2_d = float(input('What is the attack target region m?   '))

pi1 = .4
pi2 =  .4
pi1_d = 0.467
pi2_d = 0.467
xi1_d =   .75
xi2_d =   .75
xi1_d_d =   .75
xi2_d_d =   .75
r1 =   3
r2 =   8
r1_d =   3
r2_d =   8

print(pi1)
print(pi2)
print(pi1_d) 
print(pi2_d)
print(xi1_d)
print(xi2_d)
print(xi1_d_d)
print(xi2_d_d)
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
#xi1_d,2 init = 0.75
#jam density = 150vpmile
#k is 85

T = float(30/3600) # seconds -> hours
t_update = float(0.01) # updating time step; secs -> hours
#xi1_d = .55
#xi2_d = .55
k_j = 150 # jam density ranges from 185-250 per mile per lane typically; they used 60 in example
k_c = k_j / 5 # critical density is approx 1/5 kj
k_1_init = 50 # depends on initial region now
k_2_init = 100 # depends on initial region now
k_default = 85 # average overall density cars per mile for both rings (should be constant) should be k >= kj/2 to enable multivaluedness
k = k_default
# also depends on initial region
v_f =  70 # miles per hour ; free flow speed 
L = 0.25 # length of each ring in miles
# constraints
#pi1 + pi2 = 1
#xi1_d, xi2_d > 0.5


# def dk1_dt_cycle_r1(k_1, t):
# 	dk1dt = -1 * (1-xi1_d) * v_f * k_c / L  * pi1_d *T * t/(3600) +  (1-xi1_d) * v_f * k_c / L  * pi2_d *T * t/(3600)
# 	return dk1dt

# k1_0 = k_c

# t = np.array(range(0,10))
# k1 = odeint(model,k1_0,t)

# print('k(t*n)', t, ' ', k1)
# print(k_j - xi1_d * (k_j-k_c))

g1 = (1-xi1_d)*v_f/L
g2 = ((1-xi1_d)*v_f*k_c)/(L*xi1_d*(k_j-k_c))
g3 = v_f*k_c / (L*(k_j-k_c))
g4 = (1-xi2_d)*v_f/L
g5 =  ((1-xi2_d)*v_f*k_c)/(L*xi2_d*(k_j-k_c))

def dk1_dt_r1(k_1,t): 
	#return -1 * (1-xi1_d) * v_f * k_1 / L * pi1_d * T * t/(3600)
	# return (-1 * g1 * k_1) * pi1_d * T * t/(3600)
	return (-1 * g1 * k_1) * t
def dk1_dt_r2(k_1,t): 
	#return -1 * (1-xi1_d) * v_f * k_c / L * pi1_d * T * t/(3600)
	# return (-1 * g1 * k_c) * pi1_d * T * t/(3600)
	return (-1 * g1 * k_c) * t
def dk1_dt_r3(k_1,t):
	#return -1 * (k_j-k_1) * ((1-xi1_d)*v_f*k_c) / (L*xi1_d*(k_j-k_1)) * pi1_d * T * t/(3600)
	return (g2*k_1 - g2*k_j) * t
def dk1_dt_r4(k_1,t): 
	#return -1 * (v_f * k_c) / (L* (k_j-k_c)) * (k_j - 2*k+ k_1) * pi1_d * T * t/(3600)
	return ((-1 * g3 * k_1) - (g3 * (k_j - 2 * k))) * t 
def dk1_dt_r5(k_1,t): 
	#return ((1-xi2_d) * v_f) / L * (2*k-k_1)  * pi2_d * T * t/(3600)
	return ((- 1 * g4 * k_1) + (2 * g4 * k)) * t
def dk1_dt_r6(k_1,t):
	#return ((1-xi2_d) * v_f * k_c) / L * pi2_d * T * t/(3600)
	return (g4 * k_c) * t
def dk1_dt_r7(k_1,t): 
	#return  ((1-xi2_d) * v_f * k_c) / (L* xi2_d * (k_j-k_c)) * (k_j-2*k+k_1) * pi2_d * T * t/(3600)
	return ((g5 * k_1) + g5 * (k_j - 2 * k)) * t
def dk1_dt_r8(k_1,t): 
	#return  (v_f*k_c)/(L*(k_j-k_c)) * (k_j-k_1) * pi2_d * T * t/(3600)
	return ((-1 * g3 * k_1) + (g3 * k_j))  * t


# def dk1_dt_r2_6(k_1,t): 
# 	return dk1_dt_r2(k_1,t) + dk1_dt_r6(k_1,t)
# def dk1_dt_r4_7(k_1,t): 
# 	return dk1_dt_r4(k_1,t) + dk1_dt_r7(k_1,t)
# def dk1_dt_r3_6(k_1,t):
# 	return dk1_dt_r3(k_1,t) + dk1_dt_r6(k_1,t)
# def dk1_dt_r3_7(k_1,t): 
# 	return dk1_dt_r3(k_1,t) + dk1_dt_r7(k_1,t)
# def dk1_dt_r3_8(k_1,t): 
# 	return dk1_dt_r3(k_1,t) + dk1_dt_r8(k_1,t)
# def dk1_dt_r4_8(k_1,t):
# 	return dk1_dt_r4(k_1,t) + dk1_dt_r8(k_1,t)

#t = np.array(range(0,10))
#k0 = 0
#k1_r2_6 = odeint(dk1_dt_r2_6,k0,t)
# k1_r4_7 = odeint(dk1_dt_r4_7,k0,t)
# k1_r3_6 = odeint(dk1_dt_r3_6,k0,t)
# k1_r3_7 = odeint(dk1_dt_r3_7,k0,t)
# k1_r3_8 = odeint(dk1_dt_r3_8,k0,t)
# k1_r4_8 = odeint(dk1_dt_r4_8,k0,t)
t_step = 0.01/3600
def update_k1(k_1,k,R1):
	def_regions(k_1, k)
	# t = np.array(range(0,2))
	k0 = k_1
	if R1 == 1:
		# t = np.linspace(0,(int)(pi1_d * T)*120,120)
		# k1_r1 = odeint(dk1_dt_r1,k0,t)
		# k_1 = k1_r1[-1]
		t = 0
		while(t < pi1_d*T):
			k_0 = dk1_dt_r1(k0,t_step) + k0
			t = t + t_step
		return k_0
	elif R1 == 2:
		# t = np.linspace(0,(int)(pi1_d * T)*120,120)
		# k1_r2 = odeint(dk1_dt_r2,k0,t)
		# k_1 = k1_r2[-1]
		t = 0
		while(t < pi1_d*T):
			k_0 = dk1_dt_r2(k0,t_step) + k0
			t = t + t_step
		return k_0
	elif R1 == 3:
		# t = np.linspace(0,(int)(pi1_d * T)*120,120)
		# k1_r3 = odeint(dk1_dt_r3,k0,t)
		# k_1 = k1_r3[-1]
		t = 0
		while(t < pi1_d*T):
			k_0 = dk1_dt_r3(k0,t_step) + k0
			t = t + t_step
		return k_0
	elif R1 == 4:
		# t = np.linspace(0,(int)(pi1_d * T)*120,120)
		# k1_r4 = odeint(dk1_dt_r4,k0,t)
		# k_1 = k1_r4[-1]
		t = 0
		while(t < pi1_d*T):
			k_0 = dk1_dt_r4(k0,t_step) + k0
			t = t + t_step
		return k_0
	elif R1 == 5:
		# t = np.linspace(0,(int)(pi2_d * T)*120,120)
		# k1_r5 = odeint(dk1_dt_r5,k0,t)
		# k_1 = k1_r5[-1]
		t = 0
		while(t < pi2_d*T):
			k_0 = dk1_dt_r5(k0,t_step) + k0
			t = t + t_step
		return k_0
	elif R1 == 6:
		# t = np.linspace(0,(int)(pi2_d * T)*120,120)
		# k1_r6 = odeint(dk1_dt_r6,k0,t)
		# k_1 = k1_r6[-1]
		t = 0
		while(t < pi2_d*T):
			k_0 = dk1_dt_r6(k0,t_step) + k0
			t = t + t_step
		return k_0
	elif R1 == 7:
		# t = np.linspace(0,(int)(pi2_d * T)*120,120)
		# k1_r7 = odeint(dk1_dt_r7,k0,t)
		# k_1 = k1_r7[-1]
		t = 0
		while(t < pi2_d*T):
			k_0 = dk1_dt_r7(k0,t_step) + k0
			t = t + t_step
		return k_0
	elif R1 == 8:
		# t = np.linspace(0,(int)(pi2_d * T)*120,120)
		# k1_r8 = odeint(dk1_dt_r8,k0,t)
		# k_1 = k1_r8[-1]
		t = 0
		while(t < pi2_d*T):
			k_0 = dk1_dt_r8(k0,t_step) + k0
			t = t + t_step
		return k_0	
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
	r_1_k = (interval ([(k_1/2), (((1-xi1_d)*k_j - (2-xi1_d)*k_c) * k_1)/(2*k_c)]))
	global r_2_k_1
	r_2_k_1 = interval([k_c,k_j-xi1_d*(k_j-k_c)])
	global r_2_k
	r_2_k = interval ([(k_1/2), (xi1_d*k_j + (1-xi1_d)*k_c + k_1)/2])
	global r_3_k_1	
	r_3_k_1 = interval([k_j-xi1_d*(k_j-k_c), k_j])
	global r_3_k
	r_3_k = interval ([(k_1/2), ((2*xi1_d-1)/(2*xi1_d) * k_j + (k_1/(2*xi1_d)))])
	global r_4_k_1
	r_4_k_1 = interval ([0, k_j])
	global r_4_k
	r_4_k = interval([max((k_j/2 -(((1-xi1_d)*k_j - (2-xi1_d)*k_c)*k_1)/(2*k_c)), ((xi1_d*k_j) + ((1-xi1_d)*k_c) + k_1)/2, (2*xi1_d - 1)/(2*xi1_d) + (k_1)/(2*xi1_d) ),(k_1+k_j)/2])
	global r_5_k_1
	r_5_k_1 = interval([0, k_j])
	global r_5_k
	r_5_k = interval ([(k_1/2), min(((k_1 + k_c)/2), ( k_1/2 + (k_c*(k_j-k_1))/(2*(1-xi2_d)*(k_j-k_c)) ))])
	global r_6_k_1
	r_6_k_1 = interval([0, k_j-(1-xi2_d)*(k_j-k_c)])
	global r_6_k
	r_6_k = interval ([((k_1+k_c)/2), ((k_j + k_1 - xi2_d*(k_j-k_c))/2) ])
	global r_7_k_1
	r_7_k_1 = interval([0,k_j])
	global r_7_k
	r_7_k = interval ([max(((k_j + k_1 - xi2_d*(k_j-k_c))/2), ((1-2*xi2_d)*k_j + k_1)/(2*(1-xi2_d))), (k_1+k_j)/2])
	global r_8_k_1
	r_8_k_1 = interval([k_j-(1-xi2_d)*(k_j-k_c), k_j])
	global r_8_k
	r_8_k = interval ([(k_1/2) - (k_c*(k_j-k_1))/(2*(1-xi2_d)*(k_j-k_c)), ((1-2*xi2_d)*k_j + k_1)/(2*(1-xi2_d))])


# def def_regions(k_1, k):
# 	global r_1_k_1
# 	r_1_k_1 = interval([0,k_c])
# 	global r_1_k
# 	for k1 in r_1_k_1
# 		r_1_k = (interval ([(k1/2), (((1-xi1_d)*k_j - (2-xi1_d)*k_c) * k1)/(2*k_c)]) | r_1_k)
# 	global r_2_k_1
# 	r_2_k_1 = interval([k_c,k_j-xi1_d*(k_j-k_c)])
# 	global r_2_k
# 	r_2_k = interval ([(k_1/2), (xi1_d*k_j + (1-xi1_d)*k_c + k_1)/2])
# 	global r_3_k_1	
# 	r_3_k_1 = interval([k_j-xi1_d*(k_j-k_c), k_j])
# 	global r_3_k
# 	r_3_k = interval ([(k_1/2), ((2*xi1_d-1)/(2*xi1_d) * k_j + (k_1/(2*xi1_d)))])
# 	global r_4_k_1
# 	r_4_k_1 = interval ([0, k_j])
# 	global r_4_k
# 	r_4_k = interval([max((k_j/2 -(((1-xi1_d)*k_j - (2-xi1_d)*k_c)*k_1)/(2*k_c)), ((xi1_d*k_j) + ((1-xi1_d)*k_c) + k_1)/2, (2*xi1_d - 1)/(2*xi1_d) + (k_1)/(2*xi1_d) ),(k_1+k_j)/2])
# 	global r_5_k_1
# 	r_5_k_1 = interval([0, k_j])
# 	global r_5_k
# 	r_5_k = interval ([(k_1/2), min(((k_1 + k_c)/2), ( k_1/2 + (k_c*(k_j-k_1))/(2*(1-xi2_d)*(k_j-k_c)) ))])
# 	global r_6_k_1
# 	r_6_k_1 = interval([0, k_j-(1-xi2_d)*(k_j-k_c)])
# 	global r_6_k
# 	r_6_k = interval ([((k_1+k_c)/2), ((k_j + k_1 - xi2_d*(k_j-k_c))/2) ])
# 	global r_7_k_1
# 	r_7_k_1 = interval([0,k_j])
# 	global r_7_k
# 	r_7_k = interval ([max(((k_j + k_1 - xi2_d*(k_j-k_c))/2), ((1-2*xi2_d)*k_j + k_1)/(2*(1-xi2_d))), (k_1+k_j)/2])
# 	global r_8_k_1
# 	r_8_k_1 = interval([k_j-(1-xi2_d)*(k_j-k_c), k_j])
# 	global r_8_k
# 	r_8_k = interval ([(k_1/2) - (k_c*(k_j-k_1))/(2*(1-xi2_d)*(k_j-k_c)), ((1-2*xi2_d)*k_j + k_1)/(2*(1-xi2_d))])


# def def_regions(k_1, k):
# 	global r_1_k_1
# 	r_1_k_1 = interval([0,k_c])
# 	global r_1_k
# 	r_1_k = interval ([(r_1_k_1.extrema[0][0]/2), (((1-xi1_d)*k_j - (2-xi1_d)*k_c) * r_1_k_1.extrema[1][0])/(2*k_c)])
# 	global r_2_k_1
# 	r_2_k_1 = interval([k_c,k_j-xi1_d*(k_j-k_c)])
# 	global r_2_k
# 	r_2_k = interval ([(r_2_k_1.extrema[0][0]/2), (xi1_d*k_j + (1-xi1_d)*k_c + r_2_k_1.extrema[1][0])/2])
# 	global r_3_k_1	
# 	r_3_k_1 = interval([k_j-xi1_d*(k_j-k_c), k_j])
# 	global r_3_k
# 	r_3_k = interval ([(r_3_k_1.extrema[0][0]/2), ((2*xi1_d-1)/(2*xi1_d) * k_j + (r_3_k_1.extrema[1][0]/(2*xi1_d)))])
# 	global r_4_k_1
# 	r_4_k_1 = interval ([0, k_j])
# 	global r_4_k
# 	r_4_k = interval([max((k_j/2 -(((1-xi1_d)*k_j - (2-xi1_d)*k_c)*r_4_k_1.extrema[0][0])/(2*k_c)), ((xi1_d*k_j) + ((1-xi1_d)*k_c) + r_4_k_1.extrema[0][0])/2, (2*xi1_d - 1)/(2*xi1_d) + (r_4_k_1.extrema[0][0])/(2*xi1_d) ),(r_4_k_1.extrema[1][0]+k_j)/2])
# 	global r_5_k_1
# 	r_5_k_1 = interval([0, k_j])
# 	global r_5_k
# 	r_5_k = interval ([(r_5_k_1.extrema[0][0]/2), min(((r_5_k_1.extrema[1][0] + k_c)/2), ( r_5_k_1.extrema[1][0]/2 + (k_c*(k_j-r_5_k_1.extrema[1][0]))/(2*(1-xi2_d)*(k_j-k_c)) ))])
# 	global r_6_k_1
# 	r_6_k_1 = interval([0, k_j-(1-xi2_d)*(k_j-k_c)])
# 	global r_6_k
# 	r_6_k = interval ([((r_6_k_1.extrema[0][0]+k_c)/2), ((k_j + r_6_k_1.extrema[1][0] - xi2_d*(k_j-k_c))/2) ])
# 	global r_7_k_1
# 	r_7_k_1 = interval([0,k_j])
# 	global r_7_k
# 	r_7_k = interval ([max(((k_j + r_7_k_1.extrema[0][0] - xi2_d*(k_j-k_c))/2), ((1-2*xi2_d)*k_j + r_7_k_1.extrema[0][0])/(2*(1-xi2_d))), (r_7_k_1.extrema[1][0]+k_j)/2])
# 	global r_8_k_1
# 	r_8_k_1 = interval([k_j-(1-xi2_d)*(k_j-k_c), k_j])
# 	global r_8_k
# 	r_8_k = interval ([(r_8_k_1.extrema[0][0]/2) - (k_c*(k_j-r_8_k_1.extrema[0][0]))/(2*(1-xi2_d)*(k_j-k_c)), ((1-2*xi2_d)*k_j + r_8_k_1.extrema[1][0])/(2*(1-xi2_d))])

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
	# print('New k_1 after cycle is: ', k_1)
	deo = (i,j)
	return deo

# # given the combination of regions, update k_1 for an entire cycle
# def update_k1_cycle(k_1, k, deo):
# 	print('phase 1')
# 	for n in range(1, (int((pi1_d * T) / t_update))):
# 		k_1 = k_1 + max(0,update_k1(k_1,k,deo[0],1)) # update for phase 1; region i
# 	print('phase 2')
# 	for n in range(1, (int((pi2_d * T) / t_update))):
# 		k_1 = k_1 + min(k_j, update_k1(k_1,k,deo[1],2)) # udpate for phase 2; region j
# 	return k_1

# Green Time Ratio Modification Attack:
#1) same green time delta pi1 delta = pi2 delta -> same effect as changing cycle time for a given fixed state
# poinc_map_3_8 = k_j * (1-np.exp((g2-g3)*pi1*T/(3600))) + k_1_init*np.exp((g2-g3)*pi1*T/(3600))

# print(g2-g3)
# def poinc_map_3_8_delta (k_1, iters):
# 	return (k_j * (1-np.exp((g2-g3)*pi1_d*T/(3600)*iters)) + k_1*np.exp((g2-g3)*pi1_d*T/(3600)*iters))


def main():
	pi1_d = 0.467
	pi2_d = 0.467
	xi1_d = 0.6
	xi2_d = 0.6
	k1_max  = k_j - (xi1_d)*(k_j-k_c) # for it to be in region 2
	k1_min = k_c
	k_1 = k_j - (xi1_d)*(k_j-k_c) - 1
	k = (xi1_d * k_j + (2-xi1_d)*k_c) / 2
	def_regions(k_1, k)
	k_values = r_2_k & r_6_k # get the values for k that satisfy the intersection of these intervals
	print('Range of valid k values: ', k_values)
	k = k_c
	k_1_prev = 0
	count = 0
	k1_trace = []
	#k = k_1/2 - (k_c*(k_j-k_1))/(2*(1-xi2_d)*(k_j-k_c))
	#k = (k_j / 2)-1
	#for k in np.arange(k_values.extrema[0][0], k_values.extrema[1][0]):
	while return_deo(k_1,k) != (3,8):
		for k in np.arange(k_c, (k_j+k_c)/2, 2):
			for k_1 in np.arange(k_c, k_j - xi1_d*(k_j-k_c), 2):
				k1_trace = []
				if (2*k < k_1):
					break
				print('kmin is: ', k_c)
				print('kmax is: ', (k_j+k_c)/2)
				print('initial k1 is set to: ', k_1)
				print('k1 max for region 2: ', k_j - (xi1_d)*(k_j-k_c))
				print('k is set to: ', k)
				#print('the initial deo is: ', return_deo(k_1,k))
				# while return_deo(k_1,k) == (2,6):
				# 	print('running')
				# 	print('k1 should not be changing')
				# run following code on that k_1, k pair
				# modify pi2_d > pi1_d such that the effect of dk1/dt of region 6 affects the system more than 2.
				# T is still assumed to be small still; pi is also assumed to be small
				pi1_d = 0.467
				pi2_d = 0.467

				print('initial deo is: ', return_deo(k_1,k))
				if(return_deo(k_1,k) != (2,6)):
					continue;
				else:
					pi1_d = 0.3
					pi2_d = 0.7
					# i = 0
					# j = 0
					# (i,j) = return_deo(k_1,k) 
					# while (return_deo(k_1,k) == (i,j)):
					# 	(i,j) = return_deo(k_1,k) 
					# 	#print('i,j: ',i,j)
					# 	k_1 = update_k1_deo(k_1,k,i,j)
					# 	#print('running')
					# 	count = count + 1
					# 	if(count > 50):
					# 		break;
					# (i,j) =	return_deo(k_1,k)
					# print('1st final i,j: ', (i,j))

					# # pi1_d = 0.46
					# # pi2_d = 0.46

					# count = 0
					# while (return_deo(k_1,k)) == (i,j):
					# 	k_1 = update_k1_deo(k_1,k,i,j) #error occurs at (3,5) not sure why it goes to 3,5 in first place
					# 	#print('k_1: ' , k_1)
					# 	#print('running') # (2,6) -> (3,6) -> (3,5) is actually sufficient to show a drop in q; but is it correct?
					# 	count = count + 1
					# 	if(count > 50):
					# 		break;
					k1_trace.append(k_1)
					(i,j) = return_deo(k_1,k) 
					while (return_deo(k_1,k)) == (i,j):
						k_1 = update_k1_deo(k_1,k,i,j) #error occurs at (3,5) not sure why it goes to 3,5 in first place
						(i,j) = return_deo(k_1,k)
						print('2nd deo is: ', return_deo(k_1,k))
						if(return_deo(k_1,k) != (3,6)):
							print('end of 3,6: ', k_1)
							print('Not (3,6)')
							break; 
						else:
							print('start of 3,6: ', k_1)
							# pi1_d = 0.467
							# pi2_d = 0.467
						k1_trace.append(k_1)
						#print('k_1: ' , k_1)
						#print('running') # (2,6) -> (3,6) -> (3,5) is actually sufficient to show a drop in q; but is it correct?
						if(k_1_prev == k_1):
							break;
						k_1_prev = k_1

					(i,j) =	return_deo(k_1,k)
					print('2nd final i,j: ', (i,j))
					if(i,j) == (3,8):
						print('start of 3,8: ' , k_1)
						k1_trace.append(k_1)
						print('We found it!')
						break;
			(i,j) =	return_deo(k_1,k)
			if(i,j) == (3,8):
				print('We found it!')
				break;

		pi1_d = 0.467
		pi2_d = 0.467
		(i,j) = return_deo(k_1,k) 
		while (return_deo(k_1,k)) == (i,j):
			k_1 = update_k1_deo(k_1,k,i,j) #error occurs at (3,5) not sure why it goes to 3,5 in first place
			(i,j) = return_deo(k_1,k)
			print(k_1)
			if(k_1 == -1):
				break;
			k1_trace.append(k_1)
			#print('k_1: ' , k_1)
			#print('running') # (2,6) -> (3,6) -> (3,5) is actually sufficient to show a drop in q; but is it correct?
			if(k_1_prev == k_1):
				break;
			k_1_prev = k_1

		print(return_deo(k_1,k))
		# if(return_deo(k_1,k)) == (3,8):
		print(k1_trace)
		x = np.array(range(0,len(k1_trace)))
		plt.plot(x, k1_trace)  
		plt.show()
		break;

	# def_regions(k_j - (1-xi2_d)*(k_j-k_c), k)
	# # # while (return_deo(k_1,k)) == (3,5):
	# # # 	(i,j) = return_deo(k_1,k) 
	# # # 	print(i,j)
	# # # 	k_1 = update_k1_deo(k_1,k,i,j)
	# # # 	print('running')
	# # def_regions(k_1,k)
	# # print(return_deo(k_1,k))
	print('end of attack, k_1: ' , k_1)
	print('k: ', k)
	print('r1k1: ', r_1_k_1)
	print('r2k1: ', r_2_k_1)
	print('r3k1: ', r_3_k_1)
	print('r4k1: ', r_4_k_1)
	print('r5k1: ', r_5_k_1)
	print('r6k1: ', r_6_k_1)
	print('r7k1: ', r_7_k_1)
	print('r8k1: ', r_8_k_1)
	print('r1k: ', r_1_k)
	print('r2k: ', r_2_k)
	print('r3k: ', r_3_k)
	print('r4k: ', r_4_k)
	print('r5k: ', r_5_k)
	print('r6k: ', r_6_k)
	print('r7k: ', r_7_k)
	print('r8k: ', r_8_k)
	# print(k1_trace)
	# x = np.array(range(0,len(k1_trace)))
	# plt.plot(x, k1_trace)  
	# plt.show()

	#beautiful_plot(x,k1_trace)

	# def_regions(k_1,k)
	# print('the current deo after attack is: ', return_deo(k_1,k))

	# pi1_d = .46
	# pi2_d = .46
	# while return_deo(k_1,k) == (3,6):
	# 	k_1 = update_k1_deo(k_1,k,3,6)
	# 	print('running')

	# # print('the current deo after changing back is: ', return_deo(k_1,k))
	# print('k1 is set to: ', k_1)
	# print('k1 max for region 2: ', k_j - (xi1_d)*(k_j-k_c))
	# print('k is set to: ', k)
	# def_regions(k_1,k)
	# print(r_1_k_1)
	# print(r_1_k)

	# for k_1 in np.arange(1.0, k_j):
	# 	r_8_k_test = interval ([(k_1/2) + (k_c*(k_j-k_1))/(2*(1-xi2_d)*(k_j-k_c)), ((1-2*xi2_d)*k_j + k_1)/(2*(1-xi2_d))])
	# 	print('k1, r_8_k: ', k_1, r_8_k_test)

	#main()