import numpy as np
import matplotlib.pyplot as plt

#define dimensions
n = 10

def randompoint(n: int):
    '''create a random coordinate from -1 to 1 on a n-D grid'''
    p = np.random.uniform(-1.0, 1.0, n)
    return p

i = 0
in_sphere = 0
while i <= 1000000:
    #create random point
    point = randompoint(n)
    #determine if point is in the sphere
    magnitude = np.sqrt(np.sum(point**2))
    if magnitude < 1:
        in_sphere += 1
    i += 1
V = (in_sphere/i)*2**n
V_theoretical = np.pi**5/120*1
relative_error = abs(V_theoretical - V)/V_theoretical
print('The volume of a',n,'Dimensional Hypersphere is:', V)
print('The relative error is:', relative_error)

