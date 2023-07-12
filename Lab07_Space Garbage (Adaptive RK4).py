'''Code Rewritten by Yinshi Liu'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import time
#function from sample code
def rhs(r):
    """ The right-hand-side of the equations
    INPUT:
    r = [x, vx, y, vy] are floats (not arrays)
    note: no explicit dependence on time
    OUTPUT:
    1x2 numpy array, rhs[0] is for x, rhs[1] is for vx, etc"""
    M = 10.
    L = 2.

    x = r[0]
    vx = r[1]
    y = r[2]
    vy = r[3]

    r2 = x**2 + y**2
    Fx, Fy = - M * np.array([x, y], float) / (r2 * np.sqrt(r2 + .25*L**2))
    return np.array([vx, Fx, vy, Fy], float)
#initial conditions
a = 0.0
b = 10.0
h_initial = 0.01
delta = 1*10**-6
t = 0
xpoints_1 = []
vxpoints_1 = []  
ypoints_1 = []
vypoints_1 = []
time_array = []
time_step = []
r_1 = np.array([1, 0, 0, 1], float)
#start time
start_time_adaptive = time.time()
#adaptive time steps, RK4 from sample code
while t < 10:
    xpoints_1.append(r_1[0])
    vxpoints_1.append(r_1[1])
    ypoints_1.append(r_1[2])
    vypoints_1.append(r_1[3])
    dt = h_initial
    #stores the values of the last step in case of recalculation
    storage = np.copy(r_1)
    r_2 = np.copy(r_1)
    
    #RK4 for t+h
    k1 = h_initial*rhs(r_1)  
    k2 = h_initial*rhs(r_1 + 0.5*k1)
    k3 = h_initial*rhs(r_1 + 0.5*k2)
    k4 = h_initial*rhs(r_1 + k3)
    r_1 += (k1 + 2*k2 + 2*k3 + k4)/6
    #RK4 for t+h+h
    k1 = h_initial*rhs(r_1)  
    k2 = h_initial*rhs(r_1 + 0.5*k1)
    k3 = h_initial*rhs(r_1 + 0.5*k2)
    k4 = h_initial*rhs(r_1 + k3)
    r_1 += (k1 + 2*k2 + 2*k3 + k4)/6
    
    #RK4 for t+2h
    h_double = 2*h_initial
    k1_2 = h_double*rhs(r_2)  
    k2_2 = h_double*rhs(r_2 + 0.5*k1_2)  
    k3_2 = h_double*rhs(r_2 + 0.5*k2_2)
    k4_2 = h_double*rhs(r_2 + k3_2)
    r_2 += (k1_2 + 2*k2_2 + 2*k3_2 + k4_2)/6
    #calculate error
    error = 1/30*(r_2 - r_1)
    x_error = error[0]
    y_error = error[2]
    #calculate ratio rho
    p = h_initial*delta/np.sqrt(x_error**2+y_error**2)
    
    if p >= 1:
        #increase step size, keep result
        h_initial = h_initial*(min(p**0.25 ,2))
        #h_initial = h_initial*p**0.25
        t += 2*dt
        time_step.append(2*dt)
    else:
        #decrease step size to make it smaller
        h_initial = h_initial*min((p**0.25), 0.5)
        #h_initial = h_initial*p**0.25
        #recalculates RK4, overwrites result
        k1_r = h_initial*rhs(storage)  
        k2_r = h_initial*rhs(storage + 0.5*k1_r)  
        k3_r = h_initial*rhs(storage + 0.5*k2_r)
        k4_r = h_initial*rhs(storage + k3_r)
        storage += (k1_r + 2*k2_r + 2*k3_r + k4_r)/6
        r_1 = storage
        t += h_initial
        time_step.append(h_initial)
    time_array.append(t)
#end time
end_time_adaptive = time.time()
time_adaptive = end_time_adaptive - start_time_adaptive
#Regular RK4 from sample code
a = 0.0
b = 10.0
N = 10000  # let's leave it at that for now
h = (b-a)/N

tpoints = np.arange(a, b, h)
xpoints = []
vxpoints = []  # the future dx/dt
ypoints = []
vypoints = []  # the future dy/dt

# below: ordering is x, dx/dt, y, dy/dt
r = np.array([1., 0., 0., 1.], float)
#start time
start_time_regular = time.time()
for t in tpoints:
    xpoints.append(r[0])
    vxpoints.append(r[1])
    ypoints.append(r[2])
    vypoints.append(r[3])
    k1 = h*rhs(r)  # all the k's are vectors
    k2 = h*rhs(r + 0.5*k1)  # note: no explicit dependence on time of the RHSs
    k3 = h*rhs(r + 0.5*k2)
    k4 = h*rhs(r + k3)
    r += (k1 + 2*k2 + 2*k3 + k4)/6
#end time
end_time_regular = time.time()
time_regular = end_time_regular - start_time_regular
#code printouts
print('time_adaptive =', time_adaptive, 's')
print('time_regular =', time_regular, 's')
#plot the results
plt.figure(figsize = (6, 4))
plt.plot(xpoints_1, ypoints_1, '.', markersize = 3, label = 'Adaptive')
#plt.plot(xpoints, ypoints, ':', label = 'Regular, n = 10000')
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title('Trajectory using adaptive RK4 method.')
plt.axis('equal')
plt.grid()
plt.tight_layout()
plt.legend()
plt.savefig('Fig 1.1.png', dpi=500)
plt.show()

#plot time steps over time
r_sum = np.sqrt(np.array(vxpoints_1)**2 + np.array(vypoints_1)**2)
relative_r = abs(r_sum - max(r_sum))/max(r_sum)
plt.figure(figsize = (6, 4))
plt.plot(time_array, time_step, label = 'Time Step')
plt.plot(time_array, relative_r, label = 'Relative r')
plt.xlabel("Time (s)")
plt.ylabel("Step size (s)")
plt.title('Timestep size over time using adaptive RK4 method.')
plt.grid()
plt.legend()
plt.savefig('Fig 1.2.png', dpi=500)
plt.show()
