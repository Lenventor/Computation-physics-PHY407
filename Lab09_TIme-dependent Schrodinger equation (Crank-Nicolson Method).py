'''Written by Yinshi Liu'''
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
#define constants
L = 1*10**(-8)
m = 9.109*10**(-31)
sigma = L/25
kappa = 500/L
x_0 = L/5
P = 1024
dt = 1*10**(-18)
N = 3000
h_bar = 1.05*10**(-34)

p = np.arange(1, P)
a = L/P
x = p*a - L/2
psi = np.zeros(P-1, dtype = "complex_")

#normalization constant
psi_0 = 1/((2*np.pi*sigma**2)**(1/4))
#initial psi 
for i in range(1, P-1):
    psi[i] = psi_0 * np.exp(-(x[i]-x_0)**2/(4*sigma**2) + 1j*kappa*x[i])
print(psi)
plt.figure(figsize = (6,4))
plt.plot(x, psi)
plt.xlabel('x(m)')
plt.ylabel('real ψ')
plt.title('wavefunction ψ at T = 0')
plt.tight_layout()
plt.savefig('Fig 1.1.png')
#define the Hamiltonian
A = -h_bar**2/(2*m*a**2)
#V(x) = 0 in potential well
B = -2*A
vec_diag = B*np.ones(P-1)
D = np.diag(vec_diag, k=0)
sup = A*np.eye(P-1, k = 1)
sub = A*np.eye(P-1, k = -1)
H = D + sup + sub

#define time independent matrix L, R
Left = np.eye(P-1) + (dt/(2*h_bar))*1j*H
Right = np.eye(P-1) - (dt/(2*h_bar))*1j*H

#define position and time arrays
time = []
x_exp = []
prob = []
psi_total = []
prob = []
pos = []
energy = []
#begin time steps
i = 0
while i < 3000:
    #calculate psi(n+1)
    v = np.matmul(Right, psi)
    psi_new = np.linalg.solve(Left, v)
    psi = np.copy(psi_new)
    psi_total.append(psi_new)
    if i == 750:
        plt.figure(figsize = (6,4))
        plt.plot(x, psi)
        plt.xlabel('x(m)')
        plt.ylabel('real ψ')
        plt.title('wavefunction ψ at T = T/4')
        plt.tight_layout()
        plt.savefig('Fig 1.2.png')
    if i == 1500:
        plt.figure(figsize = (6,4))
        plt.plot(x, psi)
        plt.xlabel('x(m)')
        plt.ylabel('real ψ')
        plt.title('wavefunction ψ at T = T/2')
        plt.tight_layout()
        plt.savefig('Fig 1.3.png')
    if i == 2250:
        plt.figure(figsize = (6,4))
        plt.plot(x, psi)
        plt.xlabel('x(m)')
        plt.ylabel('real ψ')
        plt.title('wavefunction ψ at T = 3T/4')
        plt.tight_layout()
        plt.savefig('Fig 1.4.png')
    time.append(i*dt)
    i += 1
plt.figure(figsize = (6,4))
plt.plot(x, psi)
plt.xlabel('x(m)')
plt.ylabel('real ψ')
plt.title('wavefunction ψ at T = T')
plt.tight_layout()
plt.savefig('Fig 1.5.png')

#verify normalization
#integrate psi*conj(psi) using trap. rule
for i in range(len(psi_total)):
    psi_t = psi_total[i]
    psi_magnitude = psi_t*np.conj(psi_t)
    prob.append(np.sum(psi_magnitude)*a)
#calculate expected value of x
for i in range(len(psi_total)):
    psi_t = psi_total[i]
    position = x*psi_t*np.conj(psi_t)
    pos.append(np.sum(position)*a)
#calculate total energy
for i in range(len(psi_total)):
    psi_t = psi_total[i]
    psi_magnitude = psi_t*np.conj(psi_t)
    E = np.matmul(H, psi_magnitude)
    energy.append(np.sum(E*a))
plt.figure()
plt.plot(time, prob)
plt.ylim(0.9, 1.1)
plt.xlabel('time(s)')
plt.ylabel('Total Probability')
plt.title('Normalization of ψ over time')
plt.savefig('Fig 1.7.png')

plt.figure()
plt.plot(time, pos)
plt.xlabel('time(s)')
plt.ylabel('Expected position')
plt.title('Expected position over time')
plt.savefig('Fig 1.6.png')

plt.figure()
plt.plot(time, energy)
plt.xlabel('time(s)')
plt.ylabel('Energy (J)')
plt.title('Energy conservation over time')
plt.savefig('Fig 1.8.png')
