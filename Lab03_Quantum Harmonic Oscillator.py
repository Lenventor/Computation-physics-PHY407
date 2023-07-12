# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 23:12:18 2022
Lab 3 Q3
@author: Landon Wang
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import ones,copy,cos,tan,pi,linspace

# Define H function
def H (n,x):
    if n == 0:
        return 1
    elif n == 1:
        return 2*x
    elif n > 0:
        return 2*x*H(n-1, x)-2*(n-1)*H(n-2, x)
    else:
        return 0

# Define the wave function
def psi (n,x):
    c = (1/np.sqrt(np.math.factorial(n)*np.sqrt(np.pi)*2**n))
    return c*np.exp(-x**2/2)*H(n,x)

# Define derivative of the wave function
def d_psi (n,x):
    c = (1/np.sqrt(np.math.factorial(n)*np.sqrt(np.pi)*2**n))
    return c*np.exp(-x**2/2)*(-x*H(n,x)+2*n*H(n-1,x))

# Gaussxw by Mark Newman
def gaussxw(N):

    # Initial approximation to roots of the Legendre polynomial
    a = linspace(3,4*N-1,N)/(4*N+2)
    x = cos(pi*a+1/(8*N*N*tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta>epsilon:
        p0 = ones(N,float)
        p1 = copy(x)
        for k in range(1,N):
            p0,p1 = p1,((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)

    return x,w

# (a) #
# Defien x
l_lim = -4
r_lim = 4
step = 100
x = np.arange(l_lim, r_lim, abs(r_lim-l_lim)/step)

# Plot
N = np.arange(0, 4)
plt.figure(figsize=(6,4))
for n in N:
    plt.plot(x, psi(n, x), label = 'n = %i'% n)
plt.title('Hermonic Ocillator of n = 0, 1, 2, 3')
plt.xlabel('x')
plt.ylabel('Ψ(n,x)')
plt.legend()
plt.savefig('Figure 3.1 Hermonic Ocillator of n = 0, 1, 2, 3.png')    

# (b) #
# Define x
l_lim = -10
r_lim = 10
step = 500
x = np.arange(l_lim, r_lim, abs(r_lim-l_lim)/step)

# Plot
n = 30
plt.figure(figsize=(6,4))
plt.plot(x, psi(n, x), label = 'n = %i'% n)
plt.title('Hermonic Ocillator of n = 30')
plt.xlabel('x')
plt.ylabel('Ψ(n,x)')
plt.legend()
plt.savefig('Figure 3.2 Hermonic Ocillator of n = 30.png')

# (c) #
# def <x^2> integrating function
def average_x_sq_func (n,x):
    return ((1+x**2)/(1-x**2)**2)*(x/(1-x**2))**2*abs(psi(n,x/(1-x**2)))**2


# def <p^2> integrating function
def average_p_sq_func (n,x):
    return ((1+x**2)/(1-x**2)**2)*abs(d_psi(n,x/(1-x**2)))**2

# define N and call gaussxw
N = 100
x, w = gaussxw(N)

# Loop for every n = 1, 2, 3,..., 15
n = np.arange(16)
data = np.zeros((16, 6))
for n in n:
    X = 0
    P = 0
    for k in range(N):
        X += w[k]*average_x_sq_func(n, x[k])
        P += w[k]*average_p_sq_func(n, x[k])
    E = (X+P)/2
    rms_X = np.sqrt(X)
    rms_P = np.sqrt(P)
    data[n][0] = n
    data[n][1] = X
    data[n][2] = P
    data[n][3] = E
    data[n][4] = rms_X
    data[n][5] = rms_P
print('Part (c) data:\n   n           <x^2>       <p^2>       E           rms(x^2)    rms(p^2)\n', data)


