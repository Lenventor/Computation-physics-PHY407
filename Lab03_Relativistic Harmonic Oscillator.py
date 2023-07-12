'''Written by Yinshi Liu'''
import numpy as np
import scipy as sp
from gaussxw import *
from pylab import *
import matplotlib.pyplot as plt
#define constants
c = 3*10**8
k = 12
m = 1
x_0 = 0.01
#define g(x)
#I broke the fraction into pieces to avoid confusuion
def g(x):
    upper = k*(x_0**2-x**2)*(2*m*c**2+k*(x_0**2-x**2)/2)
    lower = 2*(m*c**2+k*(x_0**2-x**2)/2)**2
    v = c*((upper/lower)**(1/2))
    return v
#where T(x) is the antiderivative of t(x).
def t(x):
    return 4/g(x)
#create empty arrays for plotting 4/g_k and 4w_k/g_k
t_8 = np.zeros(8)
t_16 = np.zeros(16)
N_8 = []
N_16 = []
#define varibles for integration
N_1 = 8
N_2 = 16
a = 0.0
b = x_0
# Calculate the sample points and weights, then map them
x_1,w_1 = gaussxwab(N_1, a, b)
x_2,w_2 = gaussxwab(N_2, a, b)

# Perform the integration for each N value
# Added calculation for 4/g_k and 4w_k/g_k in the for loop
T_1 = 0
T_2 = 0
for i in range(N_1):
    T_1 += w_1[i]*t(x_1[i])
    t_8[i] = t(x_1[i])
    N_8.append(x_1[i])
print ('period with N=8:', T_1)
for i in range(N_2):
    T_2 += w_2[i]*t(x_2[i])
    t_16[i] = t(x_2[i])
    N_16.append(x_2[i])
print ('period with N=16:',T_2)
error_1 = 2*np.pi*np.sqrt(m/k) - T_1
error_frac_1 = error_1/2*np.pi*np.sqrt(m/k)
error_2 = 2*np.pi*np.sqrt(m/k) - T_2
error_frac_2 = error_2/2*np.pi*np.sqrt(m/k)

print ('estimated fractional error (N=8) =', error_frac_1)
print ('estimated fractional error (N=16) =', error_frac_2)
print ('classical value =', 2*np.pi*np.sqrt(m/k))

#plot for part b
plt.scatter(gaussxw(N_1)[0], t_8, label = 'N=8')
plt.scatter(gaussxw(N_2)[0], t_16, label = 'N=16')
plt.legend()
plt.xlabel('Sample position x')
plt.ylabel('value of 4/g_k')
plt.title('Value of integrand 4/g_k at each sampling point')
plt.tight_layout()
show()

plt.bar(gaussxw(N_1)[0], gaussxw(N_1)[1], width =0.02,
        label = 'N=8')
plt.bar(gaussxw(N_2)[0], gaussxw(N_2)[1], width = 0.02,
        label= 'N=16')
plt.xlabel('Sample position x')
plt.ylabel('Weight w')
plt.legend()
plt.title('Weights at each sampling point')
plt.tight_layout()
show()

plt.bar(gaussxw(N_1)[0], gaussxw(N_1)[1]*t_8, width =0.02,
        label = 'N=8')
plt.bar(gaussxw(N_2)[0], gaussxw(N_2)[1]*t_16, width = 0.02,
        label= 'N=16')
plt.xlabel('Sample position x')
plt.ylabel('Weighted value of 4w_k/g_k')
plt.legend()
plt.title('Weighted value at each sampling point')
plt.tight_layout()
show()

#part d
N_200 = 200
T_200 = 0
x_200,w_200 = gaussxwab(N_200, a, b)
for i in range(N_200):
    T_200 += w_200[i]*t(x_200[i])
print ('period with N=200:',T_200)
error_200 = (2*np.pi*np.sqrt(m/k) - T_200)/2*np.pi*np.sqrt(m/k)
print ('fractional error for N = 200:', error_200)

#plot for part e
#Define constants
x_0 = 1
x_c = 8.66*10**8
#entirely arbitary choice to speed up computation
dx = 100000
#initialize the loop
T_values = []
x_values = []
#while loop to iterate through x values
while x_0 < x_c:
    N = 16
    a = 0.0
    b = x_0
    T = 0
    x, w = gaussxw(N)
    xp = 0.5*(b-a)*x + 0.5*(b+a)
    wp = 0.5*(b-a)*w
    for i in range(N):
        T += wp[i]*t(xp[i])
    T_values.append(T)
    x_values.append(x_0)
    x_0 += dx
T_values = np.array(T_values)
x_values = np.array(x_values)

#Plots and cosmetics
plot(x_values, T_values,
     label = 'Calculated value of T')
#Calculate the classical limit
classical_T = np.ones(len(x_values))
classical_T = classical_T * 2*np.pi*np.sqrt(m/k)
#Plot the relitivistic limit
def relativistic_T (x):
    return 4*x/c
plot(x_values, relativistic_T(x_values), 'r',
     label = 'relativistic limit')
plot(x_values, classical_T, 'g',
     label = 'classical limit')
xlabel('Initial Displacement x_0 (m)')
ylabel('Time (s)')
title('Periods of Relativistic Mass Spring System')
legend()
show()


