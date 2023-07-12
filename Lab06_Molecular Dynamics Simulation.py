"""
Written by Yinshi Liu
"""
import numpy as np
import matplotlib.pyplot as plt

# define acceleration function
def a (r):
    '''
    Return the acceleration due to Lennard-Jones potential 
    when multiple particles interact

    Parameters
    ----------
    r: float
    distance to other particles

    Returns
    -------
    float
        Acceleration due to Lennard-Jones potential.

    '''
    return -48*(1/r**13)+24*(1/r**7)

#define potential function
def v(r):
    return 4*((1/r**12)-(1/r**6))

# set initial conditions according to x and y for particle 1 and 2
N = 16
Lx = 4.0
Ly = 4.0
dx = Lx/np.sqrt(N)
dy = Ly/np.sqrt(N)
x_grid = np.arange(dx/2, Lx, dx)
y_grid = np.arange(dy/2, Ly, dy)
xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)
x_initial = xx_grid.flatten()
y_initial = yy_grid.flatten()
v_initial = np.zeros(16)
acc = []
dt = 0.01
steps = 1000
accx0 = []
accy0 = []
V0 = []
#Setup Initial Conditions
for i in range(16):
    #Distance from this point to all other points
    distx = x_initial - x_initial[i]
    disty = y_initial - y_initial[i]
    #Remove the point to avoid a divide by zero error
    distx = np.delete(distx, i)
    disty = np.delete(disty, i)
    r = np.sqrt(distx**2 + disty**2)
    #calculate the angle for each interaction
    theta0 = np.arctan2(disty, distx)
    #calcualte each acceleration and split it into components
    acc = a(r)
    potential_0 = v(r)
    accx = 0; accy = 0
    vrx = 0; vry = 0
    for i in range(15):
        accx += np.cos(theta0[i])*acc[i]
        accy += np.sin(theta0[i])*acc[i]
    accx0.append(accx)
    accy0.append(accy)
    V0.append(np.sum(v(r))/2)
accx0 = np.array(accx0)
accy0 = np.array(accy0)
# calculate the initial midpoint velocity according to equation 1.10
vx_mid = dt*accx0/2
vy_mid = dt*accy0/2

#Loop for each time step
rx = []
rx.append(x_initial)
ry = []
ry.append(y_initial)
vx = []
vx.append(v_initial)
vy = []
vy.append(v_initial)
KE = []
KE.append(0.5*v_initial)
PE = []
PE.append(V0)
for i in range(1,steps):
    # calculate position according to eq. 1.11
    rxi = rx[i-1] + dt*vx_mid
    ryi = ry[i-1] + dt*vy_mid
    rx.append(rxi)
    ry.append(ryi)
    accx_step = []
    accy_step = []
    # calculate acceleration and theta at new position
    j = 0
    V = []
    while j < 16:
        rx_step = np.array(rx)
        ry_step = np.array(ry)
    #Distance from this point to all other points
        distx = rx_step[i] - rx_step[i][j]
        disty = ry_step[i] - ry_step[i][j]
    #Remove the point to avoid a divide by zero error
        distx = np.delete(distx, j)
        disty = np.delete(disty, j)
        r = np.sqrt(distx**2 + disty**2)
    #calculate the angle for each interaction
        theta0 = np.arctan2(disty, distx)
    #calcualte each acceleration and split it into components
        acc = a(r)
        potential = v(r)
        accx = 0; accy = 0
        for k in range(15):
            accx += np.cos(theta0[k])*acc[k]
            accy += np.sin(theta0[k])*acc[k]
        accx_step.append(accx)
        accy_step.append(accy)
        V.append(np.sum(v(r))/2)
        j += 1
    accx_step = np.array(accx_step)
    accy_step = np.array(accy_step)
    # calculate k as eq. 1.12
    kx = dt*accx_step
    ky = dt*accy_step
    # calculate velocity as eq. 1.13
    vxi = vx_mid + kx/2
    vyi = vy_mid + ky/2
    vx.append(vxi)
    vy.append(vyi)
    vri = vxi**2 + vyi**2
    #update energy
    KE.append(vri/2)
    PE.append(np.array(V))
    # calculate new midpoint velocity as eq. 1.14
    vx_mid = vx_mid + kx
    vy_mid = vy_mid + ky
    
E = np.array(KE) + np.array(PE)
#Extract x, y values
x1 = [item[0] for item in rx]; y1 = [item[0] for item in ry]
x2 = [item[1] for item in rx]; y2 = [item[1] for item in ry]
x3 = [item[2] for item in rx]; y3 = [item[2] for item in ry]
x4 = [item[3] for item in rx]; y4 = [item[3] for item in ry]
x5 = [item[4] for item in rx]; y5 = [item[4] for item in ry]
x6 = [item[5] for item in rx]; y6 = [item[5] for item in ry]
x7 = [item[6] for item in rx]; y7 = [item[6] for item in ry]
x8 = [item[7] for item in rx]; y8 = [item[7] for item in ry]
x9 = [item[8] for item in rx]; y9 = [item[8] for item in ry]
x10= [item[9] for item in rx]; y10= [item[9] for item in ry]
x11= [item[10] for item in rx];y11=[item[10] for item in ry]
x12= [item[11] for item in rx];y12=[item[11] for item in ry]
x13= [item[12] for item in rx];y13=[item[12] for item in ry]
x14= [item[13] for item in rx];y14=[item[13] for item in ry]
x15= [item[14] for item in rx];y15=[item[14] for item in ry]
x16= [item[15] for item in rx];y16=[item[15] for item in ry]
#Extract Energy values
E1 = [item[0] for item in E] 
E2 = [item[1] for item in E]
E3 = [item[2] for item in E]
E4 = [item[3] for item in E]
E5 = [item[4] for item in E]
E6 = [item[5] for item in E]
E7 = [item[6] for item in E]
E8 = [item[7] for item in E]
E9 = [item[8] for item in E]
E10= [item[9] for item in E]
E11= [item[10] for item in E]
E12= [item[11] for item in E]
E13= [item[12] for item in E]
E14= [item[13] for item in E]
E15= [item[14] for item in E]
E16= [item[15] for item in E]
Total = (np.array(E1)+np.array(E2)+np.array(E3)+np.array(E4)+np.array(E5)+
        np.array(E6)+np.array(E7)+np.array(E8)+np.array(E9)+np.array(E10)+
        np.array(E11)+np.array(E12)+np.array(E13)+np.array(E14)+np.array(E15)+
        np.array(E16))

conservation = (max(abs(Total))-np.mean(abs(Total)))/np.mean(abs(Total))
#plot
plt.figure(figsize=(8, 6), dpi=1000)
plt.plot(x1, y1, '.', label = 'Particle 1')
plt.plot(x2, y2, '.', label = 'Particle 2')
plt.plot(x3, y3, '.', label = 'Particle 3')
plt.plot(x4, y4, '.', label = 'Particle 4')
plt.plot(x5, y5, '.', label = 'Particle 5')
plt.plot(x6, y6, '.', label = 'Particle 6')
plt.plot(x7, y7, '.', label = 'Particle 7')
plt.plot(x8, y8, '.', label = 'Particle 8')
plt.plot(x9, y9, '.',  label = 'Particle 9')
plt.plot(x12, y12, '.',  label = 'Particle 12')
plt.plot(x13, y13, '.',  label = 'Particle 13')
plt.plot(x14, y14, '.',  label = 'Particle 14')
plt.plot(x15, y15, '.',  label = 'Particle 15')
plt.plot(x16, y16, '.',  label = 'Particle 16')
plt.plot(x10, y10, '.',  label = 'Particle 10')
plt.plot(x11, y11, '.',  label = 'Particle 11')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')
plt.title('Trajectery plot for grid of particles')
plt.legend()
plt.savefig('Figure 2.1.png') 

time_array = np.arange(0, 1000, 1)
plt.figure(figsize=(6, 4))
plt.plot(time_array, E1,  label = 'Particle 1 (Corner)')
plt.plot(time_array, E2,  label = 'Particle 2 (Side)')
plt.plot(time_array, E6,  label = 'Particle 6 (Center)')
plt.plot(time_array, Total, label = 'Total')

plt.xlabel('Time')
plt.ylabel('Energy')
plt.title('Energy for each particle')
plt.legend()
plt.savefig('Figure 2.2.png')     

print('energy is conserved within', conservation*100, '%')
