'''Written by Yinshi Liu'''
import numpy as np
import scipy as sp
from pylab import *
import matplotlib.pyplot as plt
from numpy import ones,copy,cos,tan,pi,linspace
#define constants
c = 3*10**8
k = 12
m = 1
x_0 = 1
x_c = 8.66*10**7

#Gaussxw from Newman
def gaussxw(N):
    '''Code from Lab 3'''
    # Initial approximation to roots of the Legendre polynomial
    a = linspace(3,4*N-1,N)/(4*N+2)
    x = cos(pi*a+1/(8*N*N*tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta>epsilon:
        p0 = ones(N,float)
        p1 = copy(x)
        for k in range(1,N):
            p0,p1 = p1,((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)

    return x,w

def gaussxwab(N,a,b):
    x,w = gaussxw(N)
    return 0.5*(b-a)*x+0.5*(b+a),0.5*(b-a)*w

#Find the period at x_1, x_2. x_3
def g(x):
    '''From Lab 3'''
    upper = k*(x_0**2-x**2)*(2*m*c**2+k*(x_0**2-x**2)/2)
    lower = 2*(m*c**2+k*(x_0**2-x**2)/2)**2
    v = c*((upper/lower)**(1/2))
    return v
#where T(x) is the antiderivative of t(x).
def t(x):
    '''From Lab 3'''
    return 4/g(x)

#define time array
dt = 0.0001
#2^18
time_array_1 = np.arange(0, 26.2144, dt)
time_array_2 = np.arange(0, 26.2144, dt)
#2^20
time_array_3 = np.arange(0, 104.8576, dt)
#initialize position and velocity arrays
x_1 = zeros(len(time_array_1))
v_1 = zeros(len(time_array_1))
x_2 = zeros(len(time_array_2))
v_2 = zeros(len(time_array_2))
x_3 = zeros(len(time_array_3))
v_3 = zeros(len(time_array_3))

#initial conditions
x_1[0] = x_0
v_1[0] = 0
x_2[0] = x_c
v_2[0] = 0

x_3[0] = 10*x_c
v_3[0] = 0

#initialize loop
i = 0
#Euler-Cromer for each value of x_0
while i < len(time_array_1)-1:
    v_1[i+1] = v_1[i] -k/m*x_1[i]*(1-(v_1[i]**2/c**2))**(3/2)*dt
    x_1[i+1] = x_1[i] + v_1[i+1]*dt
    i += 1
    
i = 0
while i < len(time_array_2)-1:    
    v_2[i+1] = v_2[i] -k/m*x_2[i]*(1-(v_2[i]**2/c**2))**(3/2)*dt
    x_2[i+1] = x_2[i] + v_2[i+1]*dt
    i += 1
    
i = 0
while i < len(time_array_3)-1:
    v_3[i+1] = v_3[i] -k/m*x_3[i]*(1-(v_3[i]**2/c**2))**(3/2)*dt
    x_3[i+1] = x_3[i] + v_3[i+1]*dt
    i += 1

#plot the result
plt.figure(figsize=(6,4))
plt.plot(time_array_1, x_1)
plt.xlabel('Time (s)')
plt.ylabel('Position (m)')
plt.title('Relativistic Spring System with $x_0$ = 1m')
plt.savefig('Figure 1.1 Position of Spring System with $x_0$ = 1m.png')

plt.figure(figsize=(6,4))
plt.plot(time_array_2, x_2)
plt.xlabel('Time (s)')
plt.ylabel('Position (m)')
plt.title('Relativistic Spring System with $x_0$ = $x_c$')
plt.savefig('Figure 1.2 Position of Spring System with $x_0$ = $x_c$.png')

plt.figure(figsize=(6,4))
plt.plot(time_array_3, x_3)
plt.xlabel('Time (s)')
plt.ylabel('Position (m)')
plt.title('Relativistic Spring System with $x_0$ = 10$x_c$')
plt.savefig('Figure 1.3 Position of Spring System with $x_0$ = 10$x_c$.png')

#part b
#scale the Fourier coefficients
c_1 = np.fft.rfft(x_1)
c_1_s = abs(c_1)/abs(max(c_1))
c_2 = np.fft.rfft(x_2)
c_2_s = abs(c_2)/abs(max(c_2))
c_3 = np.fft.rfft(x_3)
c_3_s = abs(c_3)/abs(max(c_3))

#convert time dimensions for each value of x0
#sampling frequency
f = int(1/dt)
N_1 = len(x_1)
N_2 = len(x_2)
N_3 = len(x_3)
t_1 = N_1*dt
t_2 = N_2*dt
t_3 = N_3*dt
#convert to angular frequency
freq_1 = np.arange(N_1/2+1)*2*pi/t_1
freq_2 = np.arange(N_2/2+1)*2*pi/t_2
freq_3 = np.arange(N_3/2+1)*2*pi/t_3
#plot the result
plt.figure(figsize=(6,4))
plt.plot(freq_1, abs(c_1_s),
         label = '$x_0$ = 1')
plt.plot(freq_2, abs(c_2_s),
         label = '$x_0$ = $x_c$')
plt.plot(freq_3, abs(c_3_s),
         label = '$x_0$ = 10$x_c$')
plt.xlabel('Angular Frequency ω (rad/s)')
plt.ylabel('Scaled Amplitude')
plt.title('Scaled Amplitudes of Fourier coefficients')
plt.legend()
plt.xlim([0, 4*2*np.pi])
plt.savefig('Figure 1.4 Scaled Amplitudes of Fourier coefficients vs. ω.png')

#part c
#Calculate the period for each x0 value
#Steps the higher the more accurate
N = 200
#define more constants
a = 0.0
b_1 = x_0
b_2 = x_c
b_3 = 10*x_c
T_1 = 0
T_2 = 0
T_3 = 0
#integration for x0, xc, 10xc
x_1, w_1 = gaussxw(N)
x_2, w_2 = gaussxw(N)
x_3, w_3 = gaussxw(N)
xp_1 = 0.5*(b_1-a)*x_1 + 0.5*(b_1+a)
xp_2 = 0.5*(b_2-a)*x_2 + 0.5*(b_2+a)
xp_3 = 0.5*(b_3-a)*x_3 + 0.5*(b_3+a)
wp_1 = 0.5*(b_1-a)*w_1
wp_2 = 0.5*(b_2-a)*w_2
wp_3 = 0.5*(b_3-a)*w_3
x_0 = 1
for i in range(N):
    T_1 += wp_1[i]*t(xp_1[i])
x_0 = x_c
for i in range(N):
    T_2 += wp_2[i]*t(xp_2[i])
x_0 = 10*x_c
for i in range(N):
    T_3 += wp_3[i]*t(xp_3[i])
#calculate frequency in Hz
f_1 = 1/T_1
f_2 = 1/T_2
f_3 = 1/T_3

#plot the result
plt.figure(figsize=(6,4))
plt.plot(freq_1/(2*np.pi), abs(c_1_s),
         label = '$x_0$ = 1')
plt.plot(freq_2/(2*np.pi), abs(c_2_s),
         label = '$x_0$ = $x_c$')
plt.plot(freq_3/(2*np.pi), abs(c_3_s),
         label = '$x_0$ = 10$x_c$')
plt.vlines(f_1,ymin = 0, ymax = 1, color = 'red',
           label = '$x_0$ = 1, calculated f = 0.5524')
plt.vlines(f_2,ymin = 0, ymax = 1, color = 'red',
           label = '$x_0$ = $x_c$, calculated f = 0.4701')
plt.vlines(f_3,ymin = 0, ymax = 1, color = 'red',
           label = '$x_0$ = 10$x_c$, calculated f = 0.0857')
plt.xlabel('Frequency f (Hz)')
plt.ylabel('Scaled Amplitude')
plt.title('Scaled Amplitudes of Fourier coefficients')
plt.legend()
plt.xlim([0, 4])
plt.savefig('Figure 1.5 Scaled Amplitudes of Fourier coefficients vs. f.png')




