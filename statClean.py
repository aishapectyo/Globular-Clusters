#!/usr/bin/env python

""" Info: This progamme compares a 'background' field with a 'cluster' field
and subtracts nearest neighbours of the bkgd stars from the cluster stars. 
By: Aisha Mahmoud-Perez """

import sys
import os
import numpy as np
from pylab import *
import scipy
import math
import matplotlib.pyplot as plt

#---Read data---#
data_in = np.genfromtxt('dat_in.txt')
ra_in = data_in[:,0]
dec_in = data_in[:,1]
t1_in = data_in[:,2]
c_in = data_in[:,3]
color_in = data_in[:,4]

data_out = np.genfromtxt('dat_out.txt')
ra_out = data_out[:,0]
dec_out = data_out[:,1]
t1_out = data_out[:,2]
c_out = data_out[:,3]
color_out = data_out[:,4]

#---Find matches between inner and outer data---#
""" Box size set to (c-t1)_out - (c-t1)_in < 0.2 and (t1_out - t_in) < 0.4
If there is a match, it exits loop. """
color_in_len = len(color_in)
color_out_len = len(color_out)
count=[]
print 'length of t1_in is '+str(len(t1_in))
#t_clean=np.zeros(len(t1_in))
#c_clean=np.zeros(len(color_in))
t_clean=[]
c_clean=[]
ra_clean=[]
dec_clean=[]

for i in xrange(len(color_in)):
	closest_match_index=-1
	closest_match_val=sys.float_info.max
	for j in xrange(len(color_out)):
		if (math.isnan(color_out[j])):
			#print 'skipped one value because it has none val'
			continue;
		if abs(color_out[j] - color_in[i]) < 0.2 and abs(t1_out[j] - t1_in[i]) < 0.4:
			temp=abs(color_out[j] - color_in[i]) + abs(t1_out[j] - t1_in[i])# distance calc			
			if (temp < closest_match_val):
				closest_match_index=j
				closest_match_val=temp				
				#print 'found another smallest similarity between color_out[j], color_in[i] ,t1_out[j] , t1_inxx[i])'
				#print 'it is as close as '+ str(closest_match_val)+ 'in index '+str(j)	
	
	if (closest_match_index== -1): #we never found a nearby star
		print ' found no nearby stars'
		print ' i value is ' + str(i) +' should be added to list with len = '+str(len(t_clean))
		t_clean.append(t1_in[i])
		c_clean.append(c_in[i])
		ra_clean.append(ra_in[i])
		dec_clean.append(dec_in[i])
		#c_clean[i]= c_in[i]		

	else:
		print 'Found adj stars remove index '+str(closest_match_index)+' for i= '+str(i)+' with closest_val= '+str(closest_match_val)
		#print 'size of colors_out '+str(len(color_out))		
		color_out[closest_match_index]=float('nan')
		t1_out[closest_match_index]=float('nan')
		color_out = np.delete(color_out, closest_match_index)
		t1_out = np.delete(t1_out, closest_match_index)
		#print 'size of colors_out '+str(len(color_out))
		#if (math.isnan(t1_out[closest_match_index])):
		#	print '         is nan baby, it is nan !'	
		#print color_out[closest_match_index]

t_clean = np.array(t_clean)
c_clean = np.array(c_clean)
ra_clean = np.array(ra_clean)
dec_clean = np.array(dec_clean)
params=np.transpose(np.array((t_clean,c_clean, ra_clean,dec_clean)))
np.savetxt('cleanData_coords.txt', params, fmt='%.14f')

plt.figure(1)
plt.plot(c_clean-t_clean, t_clean,marker='+',linestyle=' ', color='darkblue')
#plt.xlim(-0.5,3)
ax=gca()
ax.set_ylim(ax.get_ylim()[::-1])
plt.xlabel('C - R (color)')
plt.ylabel('R (magnitude)')
plt.title('Cleaned Color-Magnitude Diagram of M49')
plt.show()
