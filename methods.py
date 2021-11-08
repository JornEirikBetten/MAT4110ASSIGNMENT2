import matplotlib 
import math
matplotlib.rcParams['text.usetex'] = True
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import time

def f(x): 
	return np.pi/2.0*np.sin(np.pi/2.0*x) 

# Composite midpoint integration method
def compositeMidpoint(m, d): 
	s = 0 # initialising a variable for holding the sum
	x = np.linspace(0, 1, m) # initialising integration points
	if (m==1): 
		stepsize=1
	else: 
		stepsize=1.0/float(m) 
	# looping through all possible combinations of indexes for the 
	# values of the different dimensions. 
	for i in range(m**(d)):
		index = np.zeros(d)
		# finding the indexes for all specific combinations
		for dim in range(d): 
			index_dim = math.floor(i/m**(dim))%m
			index[dim] = index_dim
		#print(index) 
		# calculating the contribution of the index values
		ds = 1
		for j in range(0, d): 
			# calculating the value of f(x) for one value
			# then multiplying it with all previous
			# iteration for the indexes
			ds *= f((index[j]+1.0/2.0)/m)
		s += ds
			
	return s/m**d
print(np.pi/2.0*np.sin(-np.pi/4)) 
print(compositeMidpoint(2, 4))

# Monte Carlo method
def monteCarlo(M, d): 
	s = 0
	for run in range(M): 
		x = np.random.uniform(size=d) 
		ds = 1
		for dim in range(d): 
			ds *= f(x[dim])
		s += ds
	s = s/M
	return s

approximations = np.zeros(91) 
m_ = np.zeros(91)
for i in range(10, 101, 1): 
	approximations[i-10] = compositeMidpoint(i, 2) 
	m_[i-10] = i

MEDIUM = 18
print("Error when m={}: {:e}".format(m_[0], approximations[0]-1))
print("Error when m={}: {:e}".format(m_[90], approximations[90]-1))
plt.plot(np.log10(m_), np.log10(abs(approximations-1)), label=r'Absolute error')
plt.title(r'$log10(\epsilon)$ against $log10(m)$', fontsize=MEDIUM)
plt.xlabel(r'$log_{10}(m)$', fontsize=MEDIUM)
plt.ylabel(r'$log_{10}(|I_{m}^{CM}-1|)$', fontsize=MEDIUM)
plt.legend(fontsize=MEDIUM)
plt.savefig("loglogmepsilon.pdf", format="pdf") 
plt.show()

max_m = 5
approx_d4 = np.zeros(max_m)
approx_d5 = np.zeros(max_m) 
approx_d6 = np.zeros(max_m) 
approx_d7 = np.zeros(max_m) 
m = np.zeros(max_m) 

for i in range(max_m): 
	m[i] = i+1
	approx_d4[i] = abs(compositeMidpoint(i+1, 4) - 1)
	approx_d5[i] = abs(compositeMidpoint(i+1, 5) - 1)
	approx_d6[i] = abs(compositeMidpoint(i+1, 6) - 1)
	approx_d7[i] = abs(compositeMidpoint(i+1, 7) - 1)
	
	
	
plt.plot(m, approx_d4, label="d=4") 
plt.plot(m, approx_d5, label="d=5") 
plt.plot(m, approx_d6, label="d=6")
plt.plot(m, approx_d7, label="d=7") 
plt.xlabel(r'm (number of datapoints)', fontsize=MEDIUM)
plt.xticks([i for i in range(max_m)])
plt.ylabel(r'$\epsilon$', fontsize=MEDIUM)
plt.title("Absolute error vs number of datapoints, d=4, 5, 6, 7", fontsize=MEDIUM)
plt.legend(fontsize=MEDIUM)
plt.savefig("errorvsdatapointsd4567.pdf", format="pdf") 
plt.show()

"""

for i in range(max_m): 
	temp_d4 = 100
	temp_d5 = 100
	temp_d6 = 100
	temp_d7 = 100
	if (approx_d4[i] < 1e-2): 
		temp_d4 = m[i]
	if (approx_d5[i] < 1e-2): 
		temp_d5 = m[i]
	if (approx_d5[i] < 1e-2): 
		temp_d6 = m[i]
	if (approx_d7[i] < 1e-2): 
		temp_d7 = m[i]
	
print("d=4: {}".format(temp_d4))
print("d=5: {}".format(temp_d5))
print("d=6: {}".format(temp_d6))
print("d=7: {}".format(temp_d7))

"""

# Finding minimum m-values for degrees 4, 5, 6, 7
# and printing it to screen

for m in range(100): 
	if (abs(compositeMidpoint(m+1, 4)-1)<1e-2): 
		print(r'The lowest integer points m^4 to get an error smaller than 1e-2:')
		print("{}".format(m+1)) 
		break
for m in range(100): 
	if (abs(compositeMidpoint(m+1, 5)-1)<1e-2): 
		print('The lowest integer points m^5 to get an error smaller than 1e-2:')
		print("{}".format(m+1)) 
		break
for m in range(100): 
	if (abs(compositeMidpoint(m+1, 6)-1)<1e-2): 
		print('The lowest integer points m^6 to get an error smaller than 1e-2:')
		print("{}".format(m+1)) 
		break
"""		
for m in range(100): 
	if (abs(compositeMidpoint(m+1, 7)-1)<1e-2): 
		print('The lowest integer points m^7 to get an error smaller than 1e-2:')
		print("{}".format(m+1)) 
		break
"""
"""
mC_4 = np.zeros(5) 
mC_5 = np.zeros(5) 
mC_6 = np.zeros(5) 
mC_7 = np.zeros(5) 

dimensions = [4, 5, 6, 7]
points = [10**i for i in range(1, 6)]
for dimension in dimensions: 
	for i in range(len(points)): 
		mC_4[i] = monteCarlo(points[i], dimension)
		mC_5[i] = monteCarlo(points[i], dimension)
		mC_6[i] = monteCarlo(points[i], dimension)
		mC_7[i] = monteCarlo(points[i], dimension)

"""
n = 1
mC_average_4 = np.zeros(n) 
mC_average_5 = np.zeros(n) 
mC_average_6 = np.zeros(n) 
mC_average_7 = np.zeros(n) 
for i in range(n): 
	mC_average_4[i] = abs(monteCarlo(10**4, 4)-1)
for i in range(n): 
	mC_average_5[i] = abs(monteCarlo(15000, 5)-1)
for i in range(n): 
	mC_average_6[i] = abs(monteCarlo(2*10**4, 6)-1)
for i in range(n): 
	mC_average_7[i] = abs(monteCarlo(4*10**4, 7)-1)

print("With points in MC") 
print("Error d=4: {:e}".format(np.average(mC_average_4)))
print("Error d=5: {:e}".format(np.average(mC_average_5)))
print("Error d=6: {:e}".format(np.average(mC_average_6)))
print("Error d=7: {:e}".format(np.average(mC_average_7)))
"""
plt.plot(np.log10(points), abs(mC_4-1), label="d=4") 
plt.plot(np.log10(points), abs(mC_5-1), label="d=5") 
plt.plot(np.log10(points), abs(mC_6-1), label="d=6")
plt.plot(np.log10(points), abs(mC_7-1), label="d=7") 
plt.xlabel(r'$log_{10}(M)$', fontsize=MEDIUM)
plt.ylabel(r'$\epsilon$', fontsize=MEDIUM)
plt.title("Monte Carlo d=(4, 5, 6, 7)", fontsize=MEDIUM)
plt.legend(fontsize=MEDIUM)
plt.savefig("MCd4567.pdf", format="pdf") 
plt.show()
"""


# Timing Monte Carlo and Composite Midpoints method

start_time = time.time()
monteCarlo(10**4, 4) 
run_time_MC4 = time.time()-start_time

start_time = time.time()
compositeMidpoint(7, 4) 
run_time_cM4 = time.time()-start_time

start_time = time.time()
monteCarlo(15000, 5) 
run_time_MC5 = time.time()-start_time

start_time = time.time()
compositeMidpoint(8, 5) 
run_time_cM5 = time.time()-start_time

start_time = time.time()
monteCarlo(2*10**4, 6) 
run_time_MC6 = time.time()-start_time

start_time = time.time()
compositeMidpoint(8, 6) 
run_time_cM6 = time.time()-start_time

start_time = time.time()
monteCarlo(4*10**4, 7) 
run_time_MC7 = time.time()-start_time

start_time = time.time()
compositeMidpoint(9, 7) 
run_time_cM7 = time.time()-start_time


print("Run times:") 
print("MC4: %.2f" % (run_time_MC4))
print("cM4: %.2f" % (run_time_cM4))
print("MC5: %.2f" % (run_time_MC5))
print("cM5: %.2f" % (run_time_cM5))
print("MC6: %.2f" % (run_time_MC6))
print("cM6: %.2f" % (run_time_cM6))
print("MC7: %.2f" % (run_time_MC7))
print("cM7: %.2f" % (run_time_cM7))
	
