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

wb = xlsw.Workbook('data.xlsx')
ws = wb.add_worksheet()
global row
row = 0
global col
col = 0

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
		for elem in y:
			ws.write(row,col,elem)
			global col
			col = col + 1
		global row
		row = row+1
		global col
		col = 0
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
	  
	# Limit the range of the plot to only where the data is.    
	# Avoid unnecessary whitespace.    
	# plt.ylim(0, 90)    
	# plt.xlim(1968, 2014)    
	  
	# Make sure your axis ticks are large enough to be easily read.    
	# # You don't want your viewers squinting to read your plot.    
	# plt.yticks(range, [str(x) + "%" for x in range], fontsize=14)    
	# plt.xticks(fontsize=14)    
	  
	# Provide tick lines across the plot to help your viewers trace along    
	# the axis ticks. Make sure that the lines are light and small so they    
	# don't obscure the primary data lines.    
	# for y in range(10, 91, 10):    
	#     plt.plot(range(1968, 2012), [y] * len(y), "--", lw=0.5, color="black", alpha=0.3)    
	  
	# Remove the tick marks; they are unnecessary with the tick lines we just plotted.    
	# plt.tick_params(axis="both", which="both", bottom="off", top="off",    
	#                 labelbottom="on", left="off", right="off", labelleft="on")    
	  
	# Now that the plot is prepared, it's time to actually plot the data!    
	# Note that I plotted the majors in order of the highest % in the final year.    
	# majors = ['Health Professions', 'Public Administration', 'Education', 'Psychology',    
	#           'Foreign Languages', 'English', 'Communications\nand Journalism',    
	#           'Art and Performance', 'Biology', 'Agriculture',    
	#           'Social Sciences and History', 'Business', 'Math and Statistics',    
	#           'Architecture', 'Physical Sciences', 'Computer Science',    
	#           'Engineering']    
	  
	# for rank, column in enumerate(majors):    
	#     # Plot each line separately with its own color, using the Tableau 20    
	#     # color set in order.    
	#     plt.plot(gender_degree_data.Year.values,    
	#             gender_degree_data[column.replace("\n", " ")].values,    
	#             lw=2.5, color=tableau20[rank])    
	  
	#     # Add a text label to the right end of every line. Most of the code below    
	#     # is adding specific offsets y position because some labels overlapped.    
	#     y_pos = gender_degree_data[column.replace("\n", " ")].values[-1] - 0.5    
	#     if column == "Foreign Languages":    
	#         y_pos += 0.5    
	#     elif column == "English":    
	#         y_pos -= 0.5    
	#     elif column == "Communications\nand Journalism":    
	#         y_pos += 0.75    
	#     elif column == "Art and Performance":    
	#         y_pos -= 0.25    
	#     elif column == "Agriculture":    
	#         y_pos += 1.25    
	#     elif column == "Social Sciences and History":    
	#         y_pos += 0.25    
	#     elif column == "Business":    
	#         y_pos -= 0.75    
	#     elif column == "Math and Statistics":    
	#         y_pos += 0.75    
	#     elif column == "Architecture":    
	#         y_pos -= 0.75    
	#     elif column == "Computer Science":    
	#         y_pos += 0.75    
	#     elif column == "Engineering":    
	#         y_pos -= 0.25    
	  
	    # # Again, make sure that all labels are large enough to be easily read    
	    # # by the viewer.    
	    # plt.text(2011.5, y_pos, column, fontsize=14, color=tableau20[rank])    
	  
	# matplotlib's title() call centers the title on the plot, but not the graph,    
	# so I used the text() call to customize where the title goes.    
	  
	# Make the title big enough so it spans the entire plot, but don't make it    
	# so big that it requires two lines to show.    
	  
	# Note that if the title is descriptive enough, it is unnecessary to include    
	# axis labels; they are self-evident, in this plot's case.    
	# plt.text(1995, 93, "Percentage of Bachelor's degrees conferred to women in the U.S.A."    
	#        ", by major (1970-2012)", fontsize=17, ha="center")    
	  
	# Always include your data source(s) and copyright notice! And for your    
	# data sources, tell your viewers exactly where the data came from,    
	# preferably with a direct link to the data. Just telling your viewers    
	# that you used data from the "U.S. Census Bureau" is completely useless:    
	# the U.S. Census Bureau provides all kinds of data, so how are your    
	# viewers supposed to know which data set you used?    
	# plt.text(1966, -8, "Data source: nces.ed.gov/programs/digest/2013menu_tables.asp"    
	#        "\nAuthor: Randy Olson (randalolson.com / @randal_olson)"    
	#        "\nNote: Some majors are missing because the historical data "    
	#        "is not available for them", fontsize=10)    
	  
	# Finally, save the figure as a PNG.    
	# You can also save it as a PDF, JPEG, etc.    
	# Just change the file extension in this call.    
	# bbox_inches="tight" removes all the extra whitespace on the edges of your plot.    
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

T = float(30) # seconds -> hours
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
v_f = 35 # miles per hour ; free flow speed 
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
	return (-1 * g1 * k_1) * pi1_d * T * t/(3600)
def dk1_dt_r2(k_1,t): 
	#return -1 * (1-xi1_d) * v_f * k_c / L * pi1_d * T * t/(3600)
	return (-1 * g1 * k_c) * pi1_d * T * t/(3600)
def dk1_dt_r3(k_1,t):
	#return -1 * (k_j-k_1) * ((1-xi1_d)*v_f*k_c) / (L*xi1_d*(k_j-k_1)) * pi1_d * T * t/(3600)
	return (g2*k_1 - g2*k_j) * pi1_d * T * t/(3600)
def dk1_dt_r4(k_1,t): 
	#return -1 * (v_f * k_c) / (L* (k_j-k_c)) * (k_j - 2*k+ k_1) * pi1_d * T * t/(3600)
	return ((-1 * g3 * k_1) - (g3 * (k_j - 2 * k))) * pi1_d * T * t/(3600)
def dk1_dt_r5(k_1,t): 
	#return ((1-xi2_d) * v_f) / L * (2*k-k_1)  * pi2_d * T * t/(3600)
	return ((- 1 * g4 * k_1) + (2 * g4 * k)) * pi2_d * T * t/(3600)
def dk1_dt_r6(k_1,t):
	#return ((1-xi2_d) * v_f * k_c) / L * pi2_d * T * t/(3600)
	return (g4 * k_c) * pi2_d * T * t/(3600)
def dk1_dt_r7(k_1,t): 
	#return  ((1-xi2_d) * v_f * k_c) / (L* xi2_d * (k_j-k_c)) * (k_j-2*k+k_1) * pi2_d * T * t/(3600)
	return ((g5 * k_1) + g5 * (k_j - 2 * k)) * pi2_d * T * t/(3600)
def dk1_dt_r8(k_1,t): 
	#return  (v_f*k_c)/(L*(k_j-k_c)) * (k_j-k_1) * pi2_d * T * t/(3600)
	return ((-1 * g3 * k_1) + (g3 * k_j)) * pi2_d * T * t/(3600)


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

#t = np.array(range(0,10))
#k0 = 0
#k1_r2_6 = odeint(dk1_dt_r2_6,k0,t)
# k1_r4_7 = odeint(dk1_dt_r4_7,k0,t)
# k1_r3_6 = odeint(dk1_dt_r3_6,k0,t)
# k1_r3_7 = odeint(dk1_dt_r3_7,k0,t)
# k1_r3_8 = odeint(dk1_dt_r3_8,k0,t)
# k1_r4_8 = odeint(dk1_dt_r4_8,k0,t)

def update_k1(k_1,k,R1):
	def_regions(k_1, k)
	t = np.array(range(0,2))
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
		k1_r3 = odeint(dk1_dt_r3,k0,t)
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
	# if R1 == 2 and R2 == 6:
	# 	k1 = update_k1(k_1,k,R1)
	# 	k1 = update_k1(k_1,k,R2)
	# 	# k1_r2_6 = odeint(dk1_dt_r2_6,k0,t)
	# 	# k_1 = k1_r2_6[-1]
	# 	return k_1
	# elif R1 == 4 and R2 == 7:
	# 	# k1_r4_7 = odeint(dk1_dt_r4_7,k0,t)
	# 	# k_1 = k1_r4_7[-1]
	# 	return k_1
	# elif R1 == 3 and R2 == 6:
	# 	# k1_r3_6 = odeint(dk1_dt_r3_6,k0,t)
	# 	# k_1 = k1_r3_6[-1]
	# 	return k_1
	# elif R1 == 3 and R2 == 7:
	# 	# k1_r3_7 = odeint(dk1_dt_r3_7,k0,t)
	# 	# k_1 = k1_r3_7[-1]
	# 	return k_1
	# elif R1 == 3 and R2 == 8:
	# 	# k1_r3_8 = odeint(dk1_dt_r3_8,k0,t)
	# 	# k_1 = k1_r3_8[-1]
	# 	return k_1
	# elif R1 == 4 and R2 == 8:
	# 	# k1_r4_8 = odeint(dk1_dt_r4_8,k0,t)
	# 	# k_1 = k1_r4_8[-1]
	# 	return k_1
	# else:
	# 	print('error in update_k1')
	# 	return -1.0

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
	r_1_k = interval ([(k_1/2), (((1-xi1_d)*k_j - (2-xi1_d)*k_c) * k_1)/(2*k_c)])
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
	xi1_d = 0.75
	xi2_d = 0.75
	k1_max  = k_j - (xi1_d)*(k_j-k_c) # for it to be in region 2
	k1_min = k_c
	k_1 = k_j - (xi1_d)*(k_j-k_c) - 1
	k = (xi1_d * k_j + (2-xi1_d)*k_c) / 2
	def_regions(k_1, k)
	k_values = r_2_k & r_6_k # get the values for k that satisfy the intersection of these intervals

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
	k_1_min = r_3_k_1.extrema[0][0]
	k_1_max = r_3_k_1.extrema[1][0]
	# print(k_1_min)
	# print(k_1_max)
	k_1_vals = frange(130,k_1_max, step)
	k_1_vals = [130]
	# pi_vals = frange(0.2, 5.4, 0.4)
	cycle_range = range(0,40)
	# for val in pi_vals:
	# 	print('pi_d', pi1_d)
	# 	pi1_d = val
	# 	graph_all(poinc_map_3_8_delta, cycle_range, k_1_vals)
	# 	#plt.show()
	# plt.show()

	# pi_vals = frange(0.2, 5.4, 0.4)
	t_vals = frange(20,50,5)
	for val in t_vals:
		print('T_d', val)
		T = val
		graph_all(poinc_map_3_8_delta, cycle_range, k_1_vals)
		#plt.show()
	plt.show()

	wb.close()
	# k = k_default
	# for k_1 in k_1_vals:
	# 	def_regions(k_1, k_default)
	# 	T_range = range(0,4*60)
	# 	graph(poinc_map_3_8_delta, T_range, k_1)
		

	# # Plot the results
	# fig = plt.figure(1,figsize=(8, 8))
	# # # density vs time
	# ax1 = fig.add_subplot(211, projection='2d')
	# ax1.plot(x,y)
	# ax1.set_xlabel('t (hours)')
	# ax1.set_ylabel('k1(t) (vpm)')
	# #ax1.set_zlabel('pi1_d')
	# #ax1.set_title(['pi1: %f, pi2: %f, xi1: %f, and xi2: %f' %(pi1_init, pi2_init, xi1, xi2)])
	# ax2 = fig.add_subplot(212, projection='2d')
	# ax2.plot(x,y)
	# ax2.set_xlabel('t (hours)')
	# ax2.set_ylabel('k1(t) (vpm)')
	# #ax2.set_zlabel('pi1_modified')

	#ax2.set_title(['pi1: %f, pi2: %f, xi1: %f, and xi2: %f ' %(pi1_init, pi2_init, xi1, xi2)])
	#plt.show()


	#2) diff green time delta, pi1 delta != pi2 delta
	# Given that initial state is (2,6) -> show that (2,6) -> (3,6) for pi and xi. Then show the faster rate for:



	# we will test across various values for k1,k according to the initial regions (i,j)

	# get interval for [k_1 min; k_1 max] for region i

	# get interval for [k min; k max] for region i

	# get interval for  [k_1 min; k_1 max] for region j

	# get interval for [k min; k max] for region j

	# get the intersect of these two intervals for region



	#k1min = float(k_j/2)
	#k1max = float(k_j)
	# k must be greater than region 8 max and lower than region 6 minimum i.e. region 3 ^ region 8
	#kmin = float(k_j/2)
	#kmax = float((xi1*k_j+(1-xi1)*k_c+k1max) /2)


	#main()