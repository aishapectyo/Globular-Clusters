#!/usr/bin/env python

"""Info: Code to convert the initial mass to the current mass, convert mass to luminosity 
and luminosity to magnitude, then plot model lines to actual data."""

import sys
import os
import numpy as np
from pylab import *
import scipy
import matplotlib.pyplot as plt

#---Functions---#
def final_mass(formation_mass):
	"""Compute final mass of cluster.  The cluster will accrete mass, yet will lose some 
	of this mass as its trajectory through the dwarf galaxy ends and stellar evolution
	dominates the mass loss. """
	sun_mass = 1.98E30
	val = ((formation_mass)/3.0)/sun_mass
	return val

def masslum_relation(m_f):
	"""Use the Mass-Luminosity Relation to calculate the luminosity of the 
	globular clusters.  For globular clusters M/L = 2."""
	M_L = 2.0
	val = m_f/M_L
	return val

def magnitude(luminosity):
	"""Compute the magnitude in R filter for all globular clusters."""
	sun_i = 4.08
	sun_r = 4.42 #Sun's R magnitude from Lick's Observatory
	sun_L = 3.846E26
	val = -2.5*np.log10(luminosity) + sun_i
	return val

def met_conversion(metallicity):
	""" Conversion the modeled metallicity, [Fe/H] 
	to the Washington filters C, T1.
	Developed by Harris & Harris 2002 using Harris 1977 and 1999
	metallicity and color values for several globular clusters.
	B_I conversion from Harris et al. 2006"""
	val1  = 1.998 + 0.748*metallicity + 0.138*(metallicity**2)
	val2 = bi = 2.16 + 0.375*metallicity
	return val2

def distance_modulus(abs_mag):
	""" Using distance modulus to find apparent magnitude. """
	distance_m49 = 14.7E6 #pc
	distance_n4696 = 38.04E6
	val = 5*log10(distance_n4696)+(abs_mag) - 5.0
	return val

#---Read Data---#
data = np.genfromtxt('bi_o_n4696.txt')
b=data[:,0]
ifil=data[:,1]
col=data[:,2]
#data_real=np.genfromtxt('cleanData_coords.txt')
#t1=data_real[:,0]
#c=data_real[:,1]
#ra=data_real[:,2]
#dec=data_real[:,3]

data_model_blue=np.genfromtxt('4696model_metallicities_blue.txt')
feh0_blue = data_model_blue[:,0]
feh_blue = data_model_blue[:,1]
m_i_blue = data_model_blue[:,2]
m_acc_blue = data_model_blue[:,3]

data_model_red=np.genfromtxt('4696model_metallicities_red.txt')
feh0_red = data_model_red[:,0]
feh_red = data_model_red[:,1]
m_i_red = data_model_red[:,2]
m_acc_red = data_model_red[:,3]

#---Main---#
mag_r_blue=[]
m_i_len_blue=len(m_i_blue)
color_blue = []
for i in xrange(m_i_len_blue):
	m_f_blue = final_mass(m_i_blue[i])
	luminosity_blue = masslum_relation(m_f_blue)
	r_convert_blue = magnitude(luminosity_blue)
	apparent_m_blue = distance_modulus(r_convert_blue)
	mag_r_blue.append(apparent_m_blue)
	get_color_blue = met_conversion(feh_blue[i])
	color_blue.append(get_color_blue)

mag_r_red=[]
m_i_len_red=len(m_i_red)
color_red = []
for i in xrange(m_i_len_red):
	m_f_red = final_mass(m_i_red[i])
	luminosity_red = masslum_relation(m_f_red)
	r_convert_red = magnitude(luminosity_red)
	apparent_m_red = distance_modulus(r_convert_red)
	mag_r_red.append(apparent_m_red)
	get_color_red = met_conversion(feh_red[i])
	color_red.append(get_color_red)
#---Plot results---#
plt.figure(1)
#plt.plot(c-t1, t1, marker='+',linestyle=' ', color='darkolivegreen')
plt.plot(color_blue[0:10], mag_r_blue[0:10], label='[Fe/H] = -1.67', color='red')
plt.plot(color_blue[10:20], mag_r_blue[10:20], label='[Fe/H] = -1.65 ', color='black')
plt.plot(color_blue[20:30], mag_r_blue[20:30], label='[Fe/H] = -1.63', color='blue')
plt.plot(color_blue[30:40], mag_r_blue[30:40], label='[Fe/H] = -1.58', color='green')
plt.plot(color_blue[40:50], mag_r_blue[40:50], label='[Fe/H] = -1.5', color='goldenrod')
plt.plot(color_blue[50:60], mag_r_blue[50:60], label='[Fe/H] = -1.37', color='mediumpurple')
plt.plot(color_blue[60:70], mag_r_blue[60:70], label='[Fe/H] = -1.19', color='olive')
plt.plot(color_blue[70:80], mag_r_blue[70:80], label='[Fe/H] = -0.96', color='silver')
plt.plot(color_blue[80:90], mag_r_blue[80:90], label='[Fe/H] = -0.71', color='c')
plt.plot(color_blue[90:100], mag_r_blue[90:100], label='[Fe/H] = -0.43', color='hotpink')
plt.plot(col, ifil, marker='+',linestyle=' ', color='black')
plt.xlim(0.0,4.0)
plt.ylim(19, 26.5)
ax=gca()
ax.set_ylim(ax.get_ylim()[::-1])
#ax=gca()
#ax.set_ylim(ax.get_ylim()[::-1])
plt.xlabel('Color')
plt.ylabel('Magnitude (I)')
plt.title('N4696 Color-Magnitude Diagram (Blue)')
#plt.legend(loc='lower right', numpoints=1)
#plt.show()

#Red data

plt.figure(2)
#plt.plot(c-t1, t1, marker='+',linestyle=' ', color='darkolivegreen')
plt.plot(color_red[0:10], mag_r_red[0:10], label='[Fe/H] = -0.270', color='red')
plt.plot(color_red[10:20], mag_r_red[10:20], label='[Fe/H] = -0.271 ', color='black')
plt.plot(color_red[20:30], mag_r_red[20:30], label='[Fe/H] = -0.271', color='blue')
plt.plot(color_red[30:40], mag_r_red[30:40], label='[Fe/H] = -0.272', color='green')
plt.plot(color_red[40:50], mag_r_red[40:50], label='[Fe/H] = -0.274', color='goldenrod')
plt.plot(color_red[50:60], mag_r_red[50:60], label='[Fe/H] = -0.279', color='mediumpurple')
plt.plot(color_red[60:70], mag_r_red[60:70], label='[Fe/H] = -0.289', color='olive')
plt.plot(color_red[70:80], mag_r_red[70:80], label='[Fe/H] = -0.309', color='silver')
plt.plot(color_red[80:90], mag_r_red[80:90], label='[Fe/H] = -0.353', color='c')
plt.plot(color_red[90:100], mag_r_red[90:100], label='[Fe/H] = -0.455', color='hotpink')
plt.plot(col, ifil, marker='+',linestyle=' ', color='black')
plt.xlim(0.0, 4)
plt.ylim(19,26.5)
ax=gca()
ax.set_ylim(ax.get_ylim()[::-1])
plt.xlabel('Color')
plt.ylabel('Magnitude (I)')
plt.title('N4696 Color-Magnitude Diagram (Red)')
#plt.legend(loc='lower right', numpoints=1)
plt.show()


