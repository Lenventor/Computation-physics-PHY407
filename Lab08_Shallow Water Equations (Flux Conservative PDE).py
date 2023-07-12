'''Written by Yinshi Liu'''
import numpy as np
import matplotlib.pyplot as plt
#define constants
L = 1
J = 50
dx = L/J
g = 9.81
H = 0.01
t1 = 0
t2 = 1
t3 = 4
dt = 0.01
A = 0.002
sigma = 0.05
mu = 0.5
#define Force functions
def F_u (u, eta):
    '''Return F(u)'''
    return 1/2 * u**2 + g*eta

def F_eta (u, eta):
    '''Return F(η)'''
    return eta*u

#initial conditions
#define x
x = np.arange(0, 1+dx, dx)
#define u
u = np.zeros(J+1, float)
#define eta(x, 0)
eta = np.empty(J+1, float)
for i in range(0, len(eta)):
    eta[i] = (H +
              A*np.exp(-(x[i]-mu)**2/sigma**2) -
              np.mean(A*np.exp(-(x-mu)**2/sigma**2)))
t = 0
#plot t = 0 case
plt.figure(figsize = (6,4))
plt.plot(x, eta, label = 'η(x), t=0s')
plt.xlabel('Position (m)')
plt.ylabel('Wave height (m)')
plt.legend()
plt.subplots_adjust(left=0.15)
plt.title('Wave height η vs.x, t = 0s')
plt.savefig('Fig 1.1.png', dpi = 500, bbox_inches='tight')
#loop over time
u_new = np.copy(u)
eta_new = np.copy(eta)
while t < 4:
    #update u
    u_new[0] = u[0]-dt/dx*(F_u(u[1],eta[1])-F_u(u[0],eta[0]))
    for j in range(1, J):
        u_new[j] = u[j]-(dt/(2*dx))*(F_u(u[j+1],eta[j+1])-
                                   F_u(u[j-1],eta[j-1]))
    u_new[J] = u[J]-dt/dx*(F_u(u[J],eta[J])-F_u(u[J-1],eta[J-1]))
    u = np.copy(u_new)

    #update eta
    eta_new[0] = eta[0]-dt/dx*(F_eta(u[1],eta[1])-F_eta(u[0],eta[0]))
    for k in range(1, J):
        eta_new[k] = eta[k]-(dt/(2*dx))*(F_eta(u[k+1],eta[k+1])-
                                       F_eta(u[k-1],eta[k-1]))
    eta_new[J] = eta[j]-dt/dx*(F_eta(u[J],eta[J])-F_eta(u[J-1],eta[J-1]))
    eta = np.copy(eta_new)
    t += dt
    #plot 1.2 in the loop
    if abs(t-t2) <= 10**-6:
        plt.figure(figsize = (6,4))
        plt.plot(x, eta, label = 'η(x), t=1s')
        plt.xlabel('Position (m)')
        plt.ylabel('Wave height (m)')
        plt.legend()
        plt.title('Wave height η vs.x, t = 1s')
        plt.savefig('Fig 1.2.png', dpi = 500, bbox_inches='tight')
#plot fig 1.3       
plt.figure(figsize = (6,4))
plt.plot(x, eta, label = 'η(x), t=4s')
plt.xlabel('Position (m)')
plt.ylabel('Wave height (m)')
plt.legend()
plt.title('Wave height η vs.x, t = 4s')
plt.savefig('Fig 1.3.png', dpi = 500, bbox_inches='tight')
