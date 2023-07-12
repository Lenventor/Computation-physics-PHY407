# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 22:16:52 2022
Lab 7 Question 2
Adopted from squarewell.py by Mark Newman
@author: Landon Wang
"""
import numpy as np
from time import time
import matplotlib.pyplot as plt
# Constants
m = 9.1094e-31     # Mass of electron
hbar = 1.0546e-34  # Planck's constant over 2*pi
e = 1.6022e-19     # Electron charge
a = 5.2918e-11     # Bohr radius
e0 = 8.8542e-12    # Vacuum permittivity
pi = 3.1415        # Pi


# Potential function
def V(x):
    return -e**2/(4*pi*e0*x)

# Schrodinger's function for RK4
def f(r,x,E,l):
    psi = r[0]
    phi = r[1]
    fpsi = phi
    fphi = psi*(l+l**2-2*m*x**2*(V(x)-E)/hbar**2)/x**2-2*phi/x
    return np.array([fpsi,fphi],float)

# Calculate the wavefunction for a particular energy
def solve(E,l):
    psi = 0.0
    phi = 1.0
    r = np.array([psi,phi],float)

    for x in np.arange(h,r_max,h):
        k1 = h*f(r,x,E,l)
        k2 = h*f(r+0.5*k1,x+0.5*h,E,l)
        k3 = h*f(r+0.5*k2,x+0.5*h,E,l)
        k4 = h*f(r+k3,x+h,E,l)
        r += (k1+2*k2+2*k3+k4)/6

    return r[0]

# Find the energy using the secant method
def secant(n, l):
    # Name string
    name_str = str(n)
    if l == 0:
        name_str += 's'
    elif l == 1:
        name_str += 'p'
    elif l == 2:
        name_str += 'f'
    
    # Initial guess
    E1 = -15*e/n**2
    E2 = -13*e/n**2
    psi2 = solve(E1,l)

    # Secant method
    while abs(E1-E2)>target:
        psi1,psi2 = psi2,solve(E2,l)
        E1,E2 = E2,E2-psi2*(E2-E1)/(psi2-psi1)

    print("E_%s ="% name_str,E2/e,"eV")
    return E2/e


'''
# Test the significance of each initial conditions
# Orignal initial condition
print('Initial Conditions: r_max = 20a, h = 0.002a, target = e/1000')
r_max = 20*a
h = 0.002*a
target = e/1000
start = time()
E_1s = secant(1, 0)
end = time()
E_2s = secant(2, 0)
E_2p = secant(2, 1)
print('Time taken for E_1s calculation:', end - start)

# Initial condition increase r_max
print('\nInitial Conditions: r_max = 100a, h = 0.002a, target = e/1000')
r_max = 100*a
h = 0.002*a
target = e/1000
start = time()
E_1s = secant(1, 0)
end = time()
E_2s = secant(2, 0)
E_2p = secant(2, 1)
print('Time taken for E_1s calculation:', end - start)

# Initial condition decrease h
print('\nInitial Conditions: r_max = 20a, h = 0.0004a, target = e/1000')
r_max = 20*a
h = 0.0004*a
target = e/1000
start = time()
E_1s = secant(1, 0)
end = time()
E_2s = secant(2, 0)
E_2p = secant(2, 1)
print('Time taken for E_1s calculation:', end - start)

# Initial condition increase r_max
print('\nInitial Conditions: r_max = 20a, h = 0.002a, target = e/100000')
r_max = 20*a
h = 0.002*a
target = e/100000
start = time()
E_1s = secant(1, 0)
end = time()
E_2s = secant(2, 0)
E_2p = secant(2, 1)
print('Time taken for E_1s calculation:', end - start)
'''

# The best intial condition that reflects the theoretical values are chossen for further use
# Eigenstate 1s
# Initial conditions
r_max = 23*a
h = 0.002*a
target = e/1000
n = 1
l = 0
E_1s = secant(n, l)

# Recover the wavefunction with RK4
psi = 0.0
phi = 1.0
r = np.array([psi,phi],float)
X = np.arange(h,r_max,h)
Psi = np.zeros(len(X))
E = E_1s*e

for i in range (0, len(X)):
    x = X[i]
    k1 = h*f(r,x,E,l)
    k2 = h*f(r+0.5*k1,x+0.5*h,E,l)
    k3 = h*f(r+0.5*k2,x+0.5*h,E,l)
    k4 = h*f(r+k3,x+h,E,l)
    r += (k1+2*k2+2*k3+k4)/6
    Psi[i] = r[0]

# Area for normalization
A = np.sum(Psi**2)*h

# Plots
plt.figure(figsize=(8, 6), dpi=500)
plt.plot(X/a,Psi**2/A)
plt.xlabel('Radial Distance (a)')
plt.ylabel('Probability')
plt.title('Probability Density Function of R(r) for E_1s State')
plt.savefig('Figure 2.1.png')

plt.figure(figsize=(8, 6), dpi=500)
plt.plot(X/a,Psi/np.sqrt(A))
plt.xlabel('Radial Distance (a)')
plt.ylabel('R(r)')
plt.title('R(r) for E_1s State')
plt.savefig('Figure 2.2.png')


# Eigenstate 2s
# Initial conditions
r_max = 32*a
n = 2
l = 0
E_2s = secant(n, l)

# Recover the wavefunction with RK4
psi = 0.0
phi = 1.0
r = np.array([psi,phi],float)
X = np.arange(h,r_max,h)
Psi = np.zeros(len(X))
E = E_2s*e

for i in range (0, len(X)):
    x = X[i]
    k1 = h*f(r,x,E,l)
    k2 = h*f(r+0.5*k1,x+0.5*h,E,l)
    k3 = h*f(r+0.5*k2,x+0.5*h,E,l)
    k4 = h*f(r+k3,x+h,E,l)
    r += (k1+2*k2+2*k3+k4)/6
    Psi[i] = r[0]

# Area for normalization
A = np.sum(Psi**2)*h

# Plots
plt.figure(figsize=(8, 6), dpi=500)
plt.plot(X/a,Psi**2/A)
plt.xlabel('Radial Distance (a)')
plt.ylabel('Probability')
plt.title('Probability Density Function of R(r) for E_2s State')
plt.savefig('Figure 2.3.png')

plt.figure(figsize=(8, 6), dpi=500)
plt.plot(X/a,Psi/np.sqrt(A))
plt.xlabel('Radial Distance (a)')
plt.ylabel('R(r)')
plt.title('R(r) for E_2s State')
plt.savefig('Figure 2.4.png')

# Eigenstate 2p
# Initial conditions
r_max = 33*a
n = 2
l = 1
E_2p = secant(n, l)

# Recover the wavefunction with RK4
psi = 0.0
phi = 1.0
r = np.array([psi,phi],float)
X = np.arange(h,r_max,h)
Psi = np.zeros(len(X))
E = E_2p*e

for i in range (0, len(X)):
    x = X[i]
    k1 = h*f(r,x,E,l)
    k2 = h*f(r+0.5*k1,x+0.5*h,E,l)
    k3 = h*f(r+0.5*k2,x+0.5*h,E,l)
    k4 = h*f(r+k3,x+h,E,l)
    r += (k1+2*k2+2*k3+k4)/6
    Psi[i] = r[0]

# Area for normalization
A = np.sum(Psi**2)*h

# Plots
plt.figure(figsize=(8, 6), dpi=500)
plt.plot(X/a,Psi**2/A)
plt.xlabel('Radial Distance (a)')
plt.ylabel('Probability')
plt.title('Probability Density Function of R(r) for E_2p State')
plt.savefig('Figure 2.5.png')

plt.figure(figsize=(8, 6), dpi=500)
plt.plot(X/a,Psi/np.sqrt(A))
plt.xlabel('Radial Distance (a)')
plt.ylabel('R(r)')
plt.title('R(r) for E_2p State')
plt.savefig('Figure 2.6.png')