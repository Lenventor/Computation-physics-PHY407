from math import sqrt,exp
import numpy as np
from random import random,randrange, seed
import matplotlib.pyplot as plt
#Modified from Newman's salesman.py
#Modified by Yinshi Liu

# Function to calculate the magnitude of a vector
def mag(x):
    return sqrt(x[0]**2+x[1]**2)

# Function to calculate the total length of the tour
def distance():
    s = 0.0
    for i in range(N):
        s += mag(r[i+1]-r[i])
    return s

#Define seed for random points
seed(100)

# Choose N city locations and calculate the initial distance
N = 15
r = np.empty([N+1,2],float)
for i in range(N):
    r[i,0] = random()
    r[i,1] = random()
r[N] = r[0]
D = distance()

#Cooling parameters
Tmax = 10.0
Tmin = 1e-3
tau = .5*1e4

#New seed for optimization paths
seed(100)
# Main loop
t = 0
T = Tmax
while T>Tmin:

    # Cooling
    t += 1
    T = Tmax*exp(-t/tau)

    # Choose two cities to swap and make sure they are distinct
    i,j = randrange(1,N),randrange(1,N)
    while i==j:
        i,j = randrange(1,N),randrange(1,N)

    # Swap them and calculate the change in distance
    oldD = D
    r[i,0],r[j,0] = r[j,0],r[i,0]
    r[i,1],r[j,1] = r[j,1],r[i,1]
    D = distance()
    deltaD = D - oldD

    # If the move is rejected, swap them back again
    if random()>exp(-deltaD/T):
        r[i,0],r[j,0] = r[j,0],r[i,0]
        r[i,1],r[j,1] = r[j,1],r[i,1]
        D = oldD
x_values = []
y_values = []

for i in range(len(r)):
    x_values.append(r[i,0])
    y_values.append(r[i,1])

plt.figure(figsize = (6,4))
plt.scatter(x_values, y_values)
plt.plot(x_values, y_values)
plt.axis('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Traveling Salesman problem with {N} points,D = {D:.5f}'
          .format(N = N, D = D))
plt.show()
