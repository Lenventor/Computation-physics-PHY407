# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 20:51:39 2022
Lab 3 Q1
@author: Landon Wang
"""
import numpy as np
import matplotlib.pyplot as plt

# (a) #
# Define function
def f (x):
    return np.exp(np.square(x)*-1)
    
# Define forward differention method
def forward_diff(x, h):
    return (f(x+h)-f(x))/h
    
# Define array of h used for test
h = np.zeros(17)
for i in range (-16, 1):
    h[i+16] = 10 ** i

# Apply differentiation
x = 0.5
analytic = -2*f(x)*x
forward_result = np.zeros((17,2))
for i in range (0, 17):
    forward_result[i][0] = forward_diff(x, h[i])
    forward_result[i][1] = abs(forward_diff(x, h[i])-analytic)

# (b) #
# plot
plt.figure(figsize=(6,4))
plt.plot(np.log10(h), np.log10(forward_result[:, 1]), 'o-')
plt.title('Log to Log Forward Difference Error Plot')
plt.xlabel('log_10(h)')
plt.ylabel('log_10(error)')
plt.savefig('Figure 1.1 Log to Log Forward Difference Error Plot.png')


# (c) #
# Define central differention method
def central_diff(x, h):
    return (f(x+h/2)-f(x-h/2))/h

# Apply differentiation
central_result = np.zeros((17,2))
for i in range (0, 17):
    central_result[i][0] = central_diff(x, h[i])
    central_result[i][1] = abs(central_diff(x, h[i])-analytic)
    
    
# plot
plt.figure(figsize=(6,4))
plt.plot(np.log10(h), np.log10(central_result[:, 1]), 'o-')
plt.title('Log to Log Central Difference Error Plot')
plt.xlabel('log_10(h)')
plt.ylabel('log_10(error)')
plt.savefig('Figure 1.2 Log to Log Central Difference Error Plot.png')