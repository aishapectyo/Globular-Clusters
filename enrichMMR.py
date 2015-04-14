#!/usr/bin/env python

import sys
import os
import numpy as np
import scipy
import matplotlib.pyplot as plt
import math
import  scipy.integrate as si

"""Info: Code that models the enrichment of globular clusters. All free parameters.
The added metallicity is allowed to change for a range of initial masses."""

#---Constants---#
sun = 1.98E30
parsec=3.085E16
z_o = 0.017*10**(-1.104) #initial metallicity.
v = 25000 #velocity in m/s.
rho = 1.8E-19 #density of the cloud. 3/4pi(mass)/radius^(3)
time = 200.0*parsec/v #size of the cloud ~ 300pc, time to cross the cloud.
time_interval=3.15E10 #given in seconds, equivalent to 1000 years. accretion period.
time_step = int(time/time_interval) # time steps to be used in loop.

#---Functions---#
def compute_mass_acc(rho, m_i, v):
	"""Computed the mass accreted by the globular cluster as it passes
	through the dwarf galaxy.  The accreted mass is determined by the
	Bondi accretion: spherical accretion into an object. It needs
	the density of the metarial, mass of object and velocity of the
	material."""
	g=6.67E-11 #gravitational constant
	cs = 10000.0
	val = (2*rho*np.pi*(g**2)*(m_i**2))/(v**2+cs**2)**(1.5) #bondi eq.
	return val

def compute_z_new(m_acc, z_o, m_i, z1):
	"""Calculate new/final metallicity after the cloud has completed
	a single pass through the dwarf's galaxy polluted gas."""
	val = ((z_o*m_i) + (z1*m_acc))/(m_i+m_acc) #determine new metallicity
	#print(val)
	return val

#---Main---#
#Create array of abundances.
array_len=10
feh1=-3.0
z1=[]
for i in xrange(array_len):
	feh1 = feh1 + 0.3
	convert = 0.017*10**(feh1)
	z1.append(convert)

#Create array of initial masses and determine accreted mass.
m0=3000 #initial mass
m_acc=[]
m_i=[]
for j in xrange(array_len):
	m0 = m0 * 2.0
	m0_solar = m0*sun
	m_initial = m0_solar
	m_i.append(m0)
	for k in xrange(time_step):
		rate_accretion_dmdt = compute_mass_acc(rho, m0_solar, v)
		m0_solar = m0_solar + rate_accretion_dmdt*time_interval
		dm = m0_solar - m_initial
		dm = dm/sun
	m_acc.append(dm)

m_len = len(m_i)
z_iterations = len(z1)
feh_new=[]
feh0=[]
mass_i=[]
mass_f=[]
mass=np.arange(10E3*sun, 10E6*sun, 10E5*sun)
for m in xrange(z_iterations):
	for n in xrange(m_len):
		z = compute_z_new(m_acc[n], z_o, m_i[n],z1[m])
		#print(m_acc[n], z_o, m_i[n], z1[m])
		#print(z)
		convert_z0 = np.log10(z_o/0.017)
		convert_znew = np.log10(z/0.017)
		feh0.append(convert_z0)
		feh_new.append(convert_znew)
		mass_i.append(mass[n])
		mass_f.append(m_acc[n]*sun)

params=np.transpose(np.array((feh0,feh_new, mass_i, mass_f)))
np.savetxt('model_metallicities_blue.txt', params, fmt='%.14f')
print('file written, model_metallicities_blue.txt')

half=0.5
plt.figure(1)
plt.plot(mass*half, feh_new[0:10])
plt.plot(mass*half, feh_new[10:20])
plt.plot(mass*half, feh_new[20:30])
plt.plot(mass*half, feh_new[30:40])
plt.plot(mass*half, feh_new[40:50])
plt.plot(mass*half, feh_new[50:60])
plt.plot(mass*half, feh_new[60:70])
plt.plot(mass*half, feh_new[70:80])
plt.plot(mass*half, feh_new[80:90])
plt.plot(mass*half, feh_new[90:100])
plt.ylabel('Metallicity')
plt.xlabel('Mass')
plt.title('Metallicity vs. Mass')

plt.show()


