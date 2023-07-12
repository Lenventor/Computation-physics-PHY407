import numpy as np
from random import random, randrange, seed
import matplotlib.pyplot as plt
#Written by Yinshi Liu

#Define function
def f(x:float, y:float):
    return x**2 - np.cos(4*np.pi*x) + (y-1)**2

#Define Cooling Prarmeters
Tmax = 25
Tmin = 1e-6
tau = 1e4
t = 0
T = Tmax

#Define x, y
x = 2
y = 2
output = f(x, y)
x_values = []
y_values = []
t_values = []
position = []
while T > Tmin:
    # Cooling
    t_values.append(t)
    t += 1
    T = Tmax*np.exp(-t/tau)
    #update lists
    x_values.append(x)
    y_values.append(y)
    pos = (x, y)
    position.append(pos)
    #randomly move the function by dx and dy
    output = f(x, y)
    dx = np.random.normal(0, 1, 1)
    dy = np.random.normal(0, 1, 1)
    new_output = f(x+dx, y+dy)
    delta = new_output-output

    #If the move is accepted, replace x and y values
    if random()<np.exp(-delta/T):
        output = np.copy(new_output)
        x = float(x+dx)
        y = float(y+dy)
    
plt.figure(figsize = (6,4))
plt.plot(t_values, x_values, '.', markersize = 3)
plt.xlabel('Time (Steps)')
plt.ylabel('x')
plt.title('x value over time')
plt.savefig('Figure 1.1.png')
plt.figure(figsize = (6,4))
plt.plot(t_values, y_values, '.', markersize = 3)
plt.xlabel('Time (Steps)')
plt.ylabel('y')
plt.title('y value over time')
plt.savefig('Figure 1.2.png')

plt.show()

print('final x value:', x_values[-1])
print('final y value:', y_values[-1])
