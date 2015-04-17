#!/usr/bin/env python

"""Info: Code to match C and R objects using the file outputs from allstar (iraf). First, it makes a cut in chi, ignoring
anything bigger than 2, and then matches C and R coordinates down to 3 arcseconds. 
By: Aisha Mahmoud-Perez"""

import sys
import os
import numpy as np
from pylab import *
import scipy
import scipy.spatial
import matplotlib.pyplot as plt
import math

#Define variable to monitor your loop
scale = 1000

#Read file and assign columns for your data. 
data_c=np.genfromtxt('N4472_C_astrom_coords.dat')
ra_c=data_c[:,0]
dec_c=data_c[:,1]
x_c=data_c[:,2]
y_c=data_c[:,3]
mag_c=data_c[:,4]
err_c=data_c[:,5]
chi_c= data_c[:,6]

data_r=np.genfromtxt('N4472_R_astrom_coords.dat')
ra_r=data_r[:,0]
dec_r=data_r[:,1]
x_r=data_r[:,2]
y_r=data_r[:,3]
mag_r=data_r[:,4]
err_r=data_r[:,5]
chi_r=data_r[:,6]

#Variable lengths
mag_c_len = len(mag_c)
mag_r_len = len(mag_r)

#Make cut by chi-squared error. 
x_c_clean=[]
y_c_clean=[]
ra_c_clean=[]
dec_c_clean=[]
mag_c_clean=[]
err_c_clean=[]
for i in xrange(mag_c_len):
	if chi_c[i] < 2:
		x_c_clean.append(x_c[i])
		y_c_clean.append(y_c[i])
		ra_c_clean.append(ra_c[i])
		dec_c_clean.append(dec_c[i])
		mag_c_clean.append(mag_c[i])
		err_c_clean.append(err_c[i])

x_r_clean=[]
y_r_clean=[]
ra_r_clean=[]
dec_r_clean=[]
mag_r_clean=[]
err_r_clean=[]
for i in xrange(mag_r_len):
	if chi_r[i] < 2:
		x_r_clean.append(x_r[i])
		y_r_clean.append(y_r[i])
		ra_r_clean.append(ra_r[i])
		dec_r_clean.append(dec_r[i])
		mag_r_clean.append(mag_r[i])
		err_r_clean.append(err_c[i])


#Match coordinates.
ra_c_clean_len = len(ra_c_clean) #new len of file.
ra_r_clean_len = len(ra_r_clean) #new len of file.
x_c_matched=[]
y_c_matched=[]
x_r_matched=[]
y_r_matched=[]
ra_c_matched=[]
dec_c_matched=[]
ra_r_matched=[]
dec_r_matched=[]
c_matched=[]

err_c_matched=[]
r_matched=[]
err_r_matched=[]
index=[]

for j in xrange(ra_c_clean_len):
	closest_match_index_c=-1
	closest_match_val=sys.float_info.max
	if (j % scale) == 0:
		print(j)
	for k in xrange(ra_r_clean_len):
		subtract_ra = abs(ra_c_clean[j]-ra_r_clean[k])
		subtract_dec = abs(dec_c_clean[j]-dec_r_clean[k])
		if subtract_ra < 0.0003 and subtract_dec < 0.0003:
			temp = subtract_ra + subtract_dec 
			if (temp < closest_match_val):
				closest_match_index_c=j
				closest_match_index_r=k
				closest_match_val=temp
	if (closest_match_index_c== -1): 
		print 'found no nearby objects'
	else:
		print 'found closest match'
		ra_c_matched.append(ra_c_clean[closest_match_index_c])
		dec_c_matched.append(dec_c_clean[closest_match_index_c])
		ra_r_matched.append(ra_r_clean[closest_match_index_r])
		dec_r_matched.append(dec_r_clean[closest_match_index_r])
		c_matched.append(mag_c_clean[closest_match_index_c])
		err_c_matched.append(err_c_clean[closest_match_index_c])
		r_matched.append(mag_r_clean[closest_match_index_r])
		err_r_matched.append(err_r_clean[closest_match_index_r])
		x_c_matched.append(x_c_clean[closest_match_index_c])
		y_c_matched.append(y_c_clean[closest_match_index_c])
		x_r_matched.append(x_r_clean[closest_match_index_r])
		y_r_matched.append(y_r_clean[closest_match_index_r])
	
params=np.transpose(np.array((x_c_matched, y_c_matched, x_r_matched, y_r_matched, ra_c_matched, dec_c_matched, ra_r_matched, dec_r_matched, c_matched, err_c_matched, r_matched, err_r_matched)))
np.savetxt('matchedList.txt', params, fmt='%.14f')
print('file written.')
