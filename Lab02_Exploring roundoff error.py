# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 23:47:55 2022
PHY407 Lab 2 Question 4
@author: Landon Wang
"""

import numpy as np
import matplotlib.pyplot as plt



# (a) #
a = 0.98
b = 1.02
n = 500
u = np.arange(a,b,(b-a)/n)

def p(u):
    return (1-u)**8

def q(u):
    return 1-8*u+28*u**2-56*u**3+70*u**4-56*u**5+28*u**6-8*u**7+u**8

plt.figure(figsize=(16,10))
plt.plot(u, p(u), 'o', label = 'P(u)') # plot p(u)
plt.plot(u, q(u), 'o', label = 'Q(u)') # plot q(u)
plt.xlabel('u')
plt.ylabel('y')
plt.title('P(u) and Q(u) near u = 1')
plt.legend()
plt.savefig('Figure 4.1 - p and q near 1.png')

# (b) #
plt.figure(figsize=(16,10))
plt.plot(u, p(u)-q(u), 'o', label = 'P(u)') # plot p(u)-q(u)
plt.xlabel('u')
plt.ylabel('y')
plt.title('P(u) - Q(u) near u = 1')
plt.legend()
plt.savefig('Figure 4.2 - p - q near 1.png')

plt.figure(figsize=(16,10))
plt.hist(p(u)-q(u), bins = 30)
plt.xlabel('u')
plt.ylabel('Frequency')
plt.title('Distribution of P(u) - Q(u) near u = 1')
plt.savefig('Figure 4.3 - p - q near 1 distribution.png')

C = 10**-16
sd_np = np.std(p(u)-q(u))
sd_est = C*np.sqrt(n)*np.sqrt(np.mean(np.square(p(u)-q(u))))
print ('Calculated Standard Deviation =', sd_np)
print ('Estimated Roundoff Error =', sd_est)

# (c) #
a = 0.98
b = 1.
n = 150
u = np.arange(a,b,(b-a)/n)
C = 1
frac_error = C*np.sqrt(np.mean(np.square(p(u)-q(u))))/(np.sqrt(n)*np.mean(p(u)-q(u)))
print ('Fractional Error =', frac_error)

a = 0.980 # Start
b = 0.984 # End
n = 1000 # Steps
u = np.arange(a,b,(b-a)/n) # x axis
y = abs(p(u)-q(u))/abs(p(u)) # y axis
d = 50 # number of points to smoothed on each side of the center
l = (b-a)/n # step width
average_y=np.zeros(n-2*d+1) # empty array (y axis) for the smoothed points
v = np.arange(a+l*d, b-l*d, l) # empty array (x axis) for the smoothed points
for i in range(0, len(average_y)):
    average_y[i] = np.mean(y[i:i+2*d]) # smooth calcultion

plt.figure(figsize=(16,10))
plt.plot(u, y, 'o', label = '|p(u)-q(u)|/|p(u)|') # plot abs(p(u)-q(u))/abs(p(u))
plt.plot(v, average_y, 'o', label = 'Smoothed') # plot abs(p(u)-q(u))/abs(p(u))

plt.xlabel('u')
plt.ylabel('Ratio')
plt.title('|p(u)-q(u)|/|p(u)| near u = 0.98')
plt.legend()
plt.savefig('Figure 4.4 - abs fractional error near 0.98.png')


# (d) #
a = 0.98
b = 1.02
n = 500
u = np.arange(a,b,(b-a)/n)

def f(u):
    return u**8/((u**4)*(u**4))

sd = np.std(f(u))
plt.figure(figsize=(16,10))
plt.plot(u, f(u)-1, 'o', label = 'F(u)-1') # plot f(u)-1

plt.xlabel('u')
plt.ylabel('y')
plt.title('F(u)-1 near u = 1')
plt.legend()
plt.savefig('Figure 4.5 - f-1 near 1.png')

sigma = C*u





