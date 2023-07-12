# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 21:16:00 2022
PHY407 Lab 2 Question 3
@author: Landon Wang
"""

import numpy as np
import scipy.constants

n = 100 # number of steps for integration

def f (z): # define the integrating function
    return z**3/((np.exp(z/(1-z))-1)*(1-z)**5)

def Integ (a,b,n): # define the integral
    l = (b-a)/n
    #s = 0.5*f(a) + 0.5*f(b) (NOTE, This is the supposed equation but it is uncalculateable due to 0/0 and infinity)
    # Take limit to f(a) and f(b) results zero, therefore we use the nextline for further calculation
    s = 0.5*0 + 0.5*0
    for k in range (1,n):
        s += f(a+k*l)
    return l*s

pi = scipy.constants.pi
k = scipy.constants.k
h = scipy.constants.h
c = scipy.constants.c

sigma_c = (2*pi*k**4/(h**3*c**2))*Integ(0, 1, n)
sigma_sp = scipy.constants.sigma
print ('Calculated Sigma =', sigma_c)
print ('Scipy Sigma =', sigma_sp)

def W (T):
    return sigma_c*T**4