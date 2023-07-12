# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 20:10:03 2022
PHY407 Lab 6 Question 1
@author: Landon Wang
"""


import numpy as np
import matplotlib.pyplot as plt


# define acceleration function
def a (x1, x2, y1, y2):
    '''
    Return the acceleration due to Lennard-Jones potential 
    when two particles interact

    Parameters
    ----------
    x1, y1 : floot
        Position of the particle of interest.
    x2, y2: floot
        Position of the particle that interactes the particle of interest.

    Returns
    -------
    floot
        Acceleration due to Lennard-Jones potential.

    '''
    r = np.sqrt((x2-x1)**2+(y2-y2)**2)
    return (24/r**7)*(2/r**6 -1)

# set initial conditions according to x and y for particle 1 and 2
initial_x1 = (4, 4.5, 2)
initial_x2 = (5.2, 5.2, 3.5)
initial_y1 = (4, 4, 3)
initial_y2 = (4, 4, 4.4)

# loop over every set of initial conditions
for n in range (0, len(initial_x1)):
    # set initial condition in the loop
    x1 = initial_x1[n]
    x2 = initial_x2[n]
    y1 = initial_y1[n]
    y2 = initial_y2[n]

    # set dt and steps
    dt = 0.01
    steps = 100

    # empty array for velocity and position
    v1x = np.zeros(steps)
    r1x = np.zeros(steps)
    v1y = np.zeros(steps)
    r1y = np.zeros(steps)
    v2x = np.zeros(steps)
    r2x = np.zeros(steps)
    v2y = np.zeros(steps)
    r2y = np.zeros(steps)

    # store the initial condition in to the array
    r1x[0] = x1
    r1y[0] = y1
    r2x[0] = x2
    r2y[0] = y2

    # calculate accelearation use defined function
    acc = a(x1, x2, y1, y2)
    
    # calculate the angle of the acceleration. 
    # this ensures the direction of force as we expected
    # theta1 for particle 1 theta2 for particle 2
    theta1 = np.arctan((y2-y1)/(x2-x1))+np.pi
    theta2 = theta1+np.pi
    
    # calculate the initial midpoint velocity according to equation 1.10
    v1x_mid = dt*acc*np.cos(theta1)/2
    v1y_mid = dt*acc*np.sin(theta1)/2
    v2x_mid = dt*acc*np.cos(theta2)/2
    v2y_mid = dt*acc*np.sin(theta2)/2

    # loop for each step
    for i in range (1, steps):
        # calculate position according to eq. 1.11
        r1x[i] = r1x[i-1] + dt*v1x_mid
        r1y[i] = r1y[i-1] + dt*v1y_mid
        r2x[i] = r2x[i-1] + dt*v2x_mid
        r2y[i] = r2y[i-1] + dt*v2y_mid
        
        # calculate acceleration and theta at new position
        acc = a(r1x[i], r2x[i], r1y[i], r2y[i])
        theta1 = np.arctan((r2y[i]-r1y[i])/(r2x[i]-r1x[i]))+np.pi
        theta2 = theta1+np.pi
        
        # calculate k as eq. 1.12
        k1x = dt*acc*np.cos(theta1)
        k1y = dt*acc*np.sin(theta1)
        k2x = dt*acc*np.cos(theta2)
        k2y = dt*acc*np.sin(theta2)
        
        # calculate velocity as eq. 1.13
        v1x[i] = v1x_mid + k1x/2
        v1y[i] = v1y_mid + k1y/2
        v2x[i] = v2x_mid + k2x/2
        v2y[i] = v2y_mid + k2y/2
        
        # calculate new midpoint velocity as eq. 1.14
        v1x_mid = v1x_mid + k1x
        v1y_mid = v1y_mid + k1y
        v2x_mid = v2x_mid + k2x
        v2y_mid = v2y_mid + k2y

    # plot
    plt.figure(figsize=(8, 6), dpi=500)
    plt.plot(r1x, r1y, '.', label = 'Particle 1')
    plt.plot(r2x, r2y, '.', label = 'Particle 2')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Trajectery plot for initial condition %i'%(n+1))
    plt.legend()
    plt.savefig('Figure 1.%i - Trajectery for initial condition %i.png' % ((n+1), (n+1)))






