#!/usr/bin/env python

"""Code that divides cleaned data in half and plots and determines tilt slope. Need dividing color from gmm.c! """
import sys
import os
import numpy as np
from pylab import *
import scipy
import scipy.spatial
import matplotlib.pyplot as plt
import math
import scipy.stats as st
import math
from array import *
import numpy
import scipy.optimize as optimization

#---Functions---#
def line(x, a, b):
	return a*x**2 + b
def get_distance(x,y,x_center,y_center):
	d = (x-x_center)**2 + (y-y_center)**2
	d = sqrt(d)
	return d
def lower_q(mag):
	""" Computes lower quartile of data. """
	nums = sort(mag)
	try:
		low_mid = ((len(nums) - 1)/4.0)
		lq = nums[low_mid]
	except TypeError: #< There were an even amount of values.
		ceil = int(math.ceil(low_mid))
		floor = int(math.floor(low_mid))
		lq = (nums[ceil] + nums[floor])/2.0
	return lq
def upper_q(mag):
	""" Computes upper quartile of data. """
	nums_u = sort(mag)
	try:
		high_mid = (len(nums_u) - 1)*0.75
		uq = nums_u[high_mid]
	except TypeError: #< There were an even amount of values
		ceil_u = int(math.ceil_u(high_mid))
		floor_u = int(math.floor_u(high_mid))
		uq = (nums_u[ceil_u] + nums_u[floor_u])/2.0
	return uq

#---Read data---#
data=np.genfromtxt('cleanData_magnitudeCut_coords.txt')
t1=data[:,0]
c=data[:,1]
ra=data[:,2]
dec=data[:,3]
col=c-t1

#---Main---#
#Divide data set in two areas. 
ra_len = len(ra)
distance = []
for i in xrange(ra_len):
	r = get_distance(ra[i]*0.9902, dec[i], 187.4445834*0.9902, 8.0005556)
	#print(r)
	distance.append(r)

distance_max = max(distance)
dividing_distance = distance_max/2.0
#print(distance_max)
#print(dividing_distance)
ra_in=[]
dec_in=[]
t1_in=[]
c_in=[]
ra_out=[]
dec_out=[]
t1_out=[]
c_out=[]
col_in=[]
col_out=[]

for j in xrange(ra_len):
	r = get_distance(ra[j]*0.9902, dec[j], 187.4445834*0.9902, 8.0005556)
	if r < dividing_distance:
		ra_in.append(ra[j])
		dec_in.append(dec[j])
		t1_in.append(t1[j])
		c_in.append(c[j])
		col_in.append(col[j])
	else:
		ra_out.append(ra[j])
		dec_out.append(dec[j])
		t1_out.append(t1[j])
		c_out.append(c[j])
		col_out.append(col[j])

t1_in = np.array(t1_in)
c_in = np.array(c_in)
t1_out = np.array(t1_out)
c_out = np.array(c_out)
col_in = np.array(col_in)
col_out = np.array(col_out)
#Create histogram
""" Using Freedman-Diaconis rule to determine bin-size.
    max - min / h, where h = q3 - q1 (upper and lower quartiles). """

iq_in = upper_q(col_in) - lower_q(col_in)
n_in = len(col_in)
h_in = 2 * iq_in * n_in**(-1./3.)
bin_size_in = int((max(col_in) - min(col_in))/h_in)
print(bin_size_in)
params=np.transpose(np.array((col_in)))
np.savetxt('bindat_in.txt', params, fmt='%.14f')

iq_out = upper_q(col_out) - lower_q(col_out)
n_out = len(col_out)
h_out = 2 * iq_out * n_out**(-1./3.)
bin_size_out = int((max(col_out) - min(col_out))/h_out)
print(bin_size_out)
params=np.transpose(np.array((col_out)))
np.savetxt('bindat_out.txt', params, fmt='%.14f')

subplot(121)
#plt.figure(1)
plt.hist(col_in, bins=bin_size_in, histtype='step', color='r')
plt.xlabel('Color')
plt.ylabel('Number of Clusters')
plt.title('Color-binned CMD (in)')

subplot(122)
plt.hist(col_out, bins=bin_size_out, histtype='step', color='r')
plt.xlabel('Color')
plt.ylabel('Number of Clusters')
plt.title('Color-binned CMD (out)')


#Make least squares fit.
#-----inner.
sub_blue_in=[]
sub_red_in=[]
r_blue_in=[]
r_red_in=[]
col_in_len = len(col_in)
for k in xrange(col_in_len):
	sub = c_in[k] - t1_in[k]
	if sub < 1.47:
		sub_blue_in.append(sub)
		r_blue_in.append(t1_in[k])
	else:
		sub_red_in.append(sub)
		r_red_in.append(t1_in[k])


sub_blue_in = np.array(sub_blue_in)
sub_red_in = np.array(sub_red_in)
r_blue_in = np.array(r_blue_in)
r_red_in = np.array(r_red_in)
blue_res_in = optimization.curve_fit(line, r_blue_in, sub_blue_in)
red_res_in = optimization.curve_fit(line, r_red_in, sub_red_in)
print('inner: blue slope and intercept: ', blue_res_in)
print('inner: red slope and intercept: ', red_res_in)
func_blue_in = line(r_blue_in, 0.000147, 1.27)
func_red_in = line(r_red_in, 0.00031, 1.61)

#-----outer.
sub_blue_out=[]
sub_red_out=[]
r_blue_out=[]
r_red_out=[]
col_out_len = len(col_out)
for i in xrange(col_out_len):
	sub = c_out[i] - t1_out[i]
	if sub < 1.602:
		sub_blue_out.append(sub)
		r_blue_out.append(t1_out[i])

	else:
		sub_red_out.append(sub)
		r_red_out.append(t1_out[i])

sub_blue_out = np.array(sub_blue_out)
sub_red_out = np.array(sub_red_out)
r_blue_out = np.array(r_blue_out)
r_red_out = np.array(r_red_out)
blue_res_out = optimization.curve_fit(line, r_blue_out, sub_blue_out)
red_res_out = optimization.curve_fit(line, r_red_out, sub_red_out)
print('outer: blue slope and intercept: ', blue_res_out)
print('outer: red slope and intercept: ', red_res_out)
func_blue_out = line(r_blue_out, 0.0002, 1.122)
func_red_out = line(r_red_out, -0.0004, 2.04)

#Plotting.
plt.figure(2)
subplot(121)
#title('Radial Color-Magnitude Diagram')
plt.plot(col_in, t1_in, marker='+',linestyle=' ', color='black')
plt.plot(func_blue_in, r_blue_in, color='b')
plt.plot(func_red_in, r_red_in,color='r')
ax=gca()
ax.set_ylim(ax.get_ylim()[::-1])
plt.xlim(0,3)
plt.xlabel('C-R (color)')
plt.ylabel('R (magnitude)')
plt.title('Color-Magnitude Diagram (in)')

subplot(122)
plt.plot(col_out, t1_out, marker='+',linestyle=' ', color='black')
plt.plot(func_blue_out, r_blue_out,color='b')
plt.plot(func_red_out, r_red_out, color='r')
ax=gca()
ax.set_ylim(ax.get_ylim()[::-1])
plt.xlim(0,3)
plt.xlabel('C-R (color)')
plt.ylabel('R (magnitude)')
plt.title('Color-Magnitude Diagram (out)')
plt.show()
