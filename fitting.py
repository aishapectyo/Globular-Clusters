import sys
import scipy
import scipy.stats
import scipy.optimize
import sys
import os
import numpy as np
from pylab import *
from array import *
import matplotlib.pyplot as plt
import math
import random
from random import randrange, uniform

"""Code that fit a function to data and computes the error. Aisha Mahmoud-Perez."""

#-------------Functions--------------#
def abs_mag(mag,d):
	"""Function that computes the absolute magnitude for the 
	globular clusters in NGC 4472 using 17.2 as the distance to the galaxy"""
	M = [mag[i] -  5*(log(d) - 1) for i in range(len(mag))]
	return M


def gauss_fit(data, nBins = 50):
	"""Fit function to data.  Computes best fit of a fit-function to a 
	pdf of the input data using the least-squares method.
	Input: data = absolute magnitudes.
	nBins = number of bins for the frequency histogram. Returns: s = resulting fit-parameters.
	"""

	#prepare observed x,y-values (bin centers and probability densities).
	N = len(data)
	#nBins = linspace(1.5,6,50)

	n__, bins__, patches__ = hist(M, nBins, color='black', label='Data')
	xVals = py.linspace(bins__[0], bins__[-1], (len(bins__)-1))
	yVals = [n__[i] for i in range(nBins)]
	#define initial guesses and fit function.
	#define objective function as the vertical difference
	#between the observed data and the fit-function
	fitFunc = lambda s,x: s[0]  * py.exp(-(x - s[1])**2.0 / (2 * s[2]**2))
	objFunc = lambda s,x,y: (fitFunc(s,x)-y)
	s0 = [15.0, 10.0, 2.0]
	s, flag = scipy.optimize.leastsq(objFunc,s0,args=(xVals,yVals))
	x_fit = py.linspace(xVals[0], xVals[-1], 49)
	y_fit = fitFunc(s, x_fit)
	plt.plot(x_fit, y_fit, lw=4, color='r', label='Gaussian fit')
	plt.xlabel('Binned Magnitudes')
	plt.ylabel('Events')
	plt.title('NGC 4472 Globular Cluster Luminosity Function')
	plt.legend(loc='upper left', numpoints = 1)
	plt.show()
	return s

def t_5(data, nBins=50):
	"""Fit t_5, Secker 1992, function to data. Computes best fit to a pdf of the input
	data using the least-squares method.
	Input: data = absolute magnitudes.
	nBins = number of bins for the frequency histogram. 
	Returns: f = resulting fit-parameters"""
	
	#prepare observed x,y-values (bin centers and probability densities).
	N = len(data)
	n__, bins__, patches__ = hist(M, nBins, color='black', label='Data')
	xVals = py.linspace(bins__[0], bins__[-1], (len(bins__)-1))
	yVals = [n__[i] for i in range(nBins)]
	
	#define initial guesses and fit function.
	fitFunc = lambda f,x: f[0] * (1 + ((x-f[1])**2/(5*f[2]**2)))**(-3)
	objFunc = lambda f,x,y: (fitFunc(f,x)-y)
	f0 = [15. ,10.0, 2.0]
	f,flag = scipy.optimize.leastsq(objFunc,f0,args=(xVals,yVals))
	x_fit = py.linspace(xVals[0], xVals[-1], nBins)
	y_fit = fitFunc(f, x_fit)
	y_init = fitFunc(f0,x_fit)
	plt.plot(x_fit, y_fit, lw=4, color='r', label = 't_5 fit')
	plt.xlabel('Binned Magnitudes')
	plt.ylabel('Events')
	plt.title('NGC 4472 Globular Cluster Luminosity Function')
	plt.legend(loc='upper left', numpoints = 1)
	plt.show()
	return f


#----------------Load data----------------#
data = py.loadtxt('matched_CR_psf.txt')
RA_C = data[:,0]
DEC_C =  data[:,1]
x_C =  data[:,2]
y_C = data[:,3]
mag_C = data[:,4]
merr_C = data[:,5]
RA_R = data[:,6]
DEC_R =  data[:,7]
x_R =  data[:,8]
y_R = data[:,9]
mag_R = data[:,10]
merr_R = data[:,11]

#-----------Call your functions-----------#
M = abs_mag(mag_C, 17.2)
r = gauss_fit(M)
mean_g = r[1]
sigma_g = r[2]
print('Gaussian values: ', mean_g, sigma_g)

r_2 = t_5(M)
mean_t = r_2[1]
sigma_t = r_2[2]
print('t_5  values: ', mean_t, sigma_t)

