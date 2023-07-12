"""
This program simulates Brownian motion in the presence of walls
Note that the physical behaviour would be to stick to walls,
which is the purpose of Q1a.
Author: Nico Grisouard, University of Toronto
"""
#Modified by Yinshi Liu
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


def nextmove(x, y):
    """ randomly choose a direction
    0 = up, 1 = down, 2 = left, 3 = right"""
    direction =  np.random.randint(4)# COMPLETE

    if direction == 0:  # move up
        y += 1
    elif direction == 1:  # move down
        y -= 1
    elif direction == 2:  # move right
        x += 1
    elif direction == 3:  # move left
        x -= 1
    else:
        print("error: direction isn't 0-3")

    return x, y



# %% main program starts here ------------------------------------------------|
# YOU NEED TO FINISH IT!

Lp = 101  # size of domain
Nt = 5000  # number of time steps
# arrays to record the trajectory of the particle
x = []
y = []
t = np.arange(0, Nt+1)
# COMPLETE

centre_point = (Lp-1)//2  # middle point of domain
xp = centre_point
yp = centre_point
x.append(xp)
y.append(yp)
i = 0
while i < Nt:
    xpp, ypp = nextmove(xp, yp)
    xp = np.copy(xpp)
    yp = np.copy(ypp)
    i += 1
    #check if particle hits wall
    while xp not in range(1, Lp-1) or yp not in range(1, Lp-1):
        #find new direction
        xpp, ypp = nextmove(xp, yp)
        i += 1
        #update position with wall limit
        xpp = max(min(xpp, Lp-1), 0)
        ypp = max(min(ypp, Lp-1), 0)
        xp = np.copy(xpp)
        yp = np.copy(ypp)
        x.append(xp)
        y.append(yp)
    x.append(xp)
    y.append(yp)
    
x = np.array(x)
y = np.array(y)
for i in range(len(x)):
    x[i] = max(min(x[i], Lp-1), 0)
    y[i] = max(min(y[i], Lp-1), 0)
#plot the walls
bottom = np.zeros(Lp)
top = bottom + Lp-1
x_box = np.arange(0, Lp)
left = np.zeros(Lp)
right = left + Lp-1
y_box = np.arange(0, Lp)

plt.figure(figsize = (6,4))
plt.plot(x, y)
plt.plot(x_box, bottom, color = 'blue')
plt.plot(x_box, top, color = 'blue')
plt.plot(top, y_box, color = 'blue')
plt.plot(bottom, y_box, color = 'blue')
plt.xlabel('i')
plt.ylabel('j')
plt.axis('equal')
plt.title('Trajectory of particle undergoing Brownian motion')
plt.savefig('Fig 1.2.png')

plt.figure(figsize = (6,4))
plt.plot(t, x)
plt.xlabel('Time (steps)')
plt.ylabel('i')
plt.title('i coordinate of particle over time')
plt.savefig('Fig 1.3.png')

plt.figure(figsize = (6,4))
plt.plot(t, y)
plt.xlabel('Time (steps)')
plt.ylabel('j')
plt.title('j coordinate of particle over time')
plt.savefig('Fig 1.4.png')
