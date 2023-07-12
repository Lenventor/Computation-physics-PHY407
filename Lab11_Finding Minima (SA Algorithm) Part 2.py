import numpy as np
from random import random, randrange, seed
import matplotlib.pyplot as plt
#Written by Yinshi Liu

#Define function
def f(x:float, y:float):
    return np.cos(x)+np.cos(np.sqrt(2)*x)+np.cos(np.sqrt(3)*x)+(y-1)**2

#Define Cooling Prarmeters
Tmax = 20
Tmin = 1e-7
tau = 1e5
t = 0
T = Tmax

#Define x, y
x = 2
y = 2
output = f(x, y)
x_values = []
y_values = []
t_values = []
while T > Tmin:
    # Cooling
    t_values.append(t)
    t += 1
    T = Tmax*np.exp(-t/tau)
    #update lists
    x_values.append(x)
    y_values.append(y)
    #randomly move the function by dx and dy
    output = f(x, y)
    dx = np.random.normal(0, 1, 1)
    dy = np.random.normal(0, 1, 1)
    new_output = f(x+dx, y+dy)
    delta = new_output-output
    x_new = float(x+dx)
    y_new = float(y+dy)
    #If the move is accepted, keep the current move
    if random()<np.exp(-delta/T):
        #check if new move is in range
        if x_new < 51 and x_new > 0 and y_new < 20 and y_new > -20:
            output = np.copy(new_output)
            x = x_new
            y = y_new
plt.figure(figsize = (6,4))
plt.plot(t_values, x_values, '.', markersize = 3)
plt.xlabel('Time (Steps)')
plt.ylabel('x')
plt.title('x value over time')
plt.savefig('Figure 1.3.png')
plt.figure(figsize = (6,4))
plt.plot(t_values, y_values, '.', markersize = 3)
plt.xlabel('Time (Steps)')
plt.ylabel('y')
plt.title('y value over time')
plt.savefig('Figure 1.4.png')
plt.show()

print('final x value:', x_values[-1])
print('final y value:', y_values[-1])
