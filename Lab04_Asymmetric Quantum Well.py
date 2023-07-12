'''Written by Yinshi Liu'''
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
#part a
#define function
def f(c,x):
    return 1-np.exp(-c*x)
#initial condition
x = 1.0
c = 2
error = 1.0
accuracy = 10**-6
#set up x
x_list = [x]
i = 0
#loop
while error > accuracy:
    x_list.append(f(c, x_list[-1]))
#Since function is known, use e = (x-x')/(1-1/f'(x))
    error = abs((x_list[-1]-x_list[-2])/(1-np.exp(c*x_list[-1])/c))
    i += 1
#print result
print ('Solution of x:', x_list[-1])
print ('Number of iterations:', i)
#part b
#new initial conditions
c_max = 3
dc = .01
nt = int(c_max/dc)
#set up arrays
x_value = []
#initialize loop
c = np.linspace(dc, c_max, nt)
for C in c:
    m1 = 1
    error = 1
    while error > accuracy:
        m1, m2 = f(C,m1), m1
        error = abs((m1 - m2)/(1-np.exp(C*m2)/C))
    x_value.append(m1)
#plot the result
plt.plot(c, x_value)
plt.xlabel('c')
plt.ylabel('x = 1-e^(-cx)')
plt.title('Solution of x from c = 0 to c = 3, dc = 0.01')
plt.show()
#set up overrelaxation
x = 1.0
c = 2
omega = 0.50
x1 = 1
error = 1
i = 0
#adapted code from part a
accuracy = 10**-6
while error > accuracy:
    dx = f(c, x1) - x1
    x1, x2 = x1+(1+omega)*dx, x1
    error = abs((x1 - x2)/(1-1/((1+omega)*c*np.exp(-c*x2)-omega)))
    i += 1
#print outputs
print ('Solution of x:', x1)
print ('Number of iterations:', i)

#Part c
#Define function
def g(x):
    return 5*np.exp(-x)+x-5
#define x1, x2 and accuracy
x_1 = 4
x_2 = 5
accuracy = 10**-6
error = abs(x_2 - x_1)
#check signs
sign_x_1 = np.sign(g(x_1))
sign_x_2 = np.sign(g(x_2))
if sign_x_1 == sign_x_2:
    print('Same signs encountered!')
#binary search
while error > accuracy:
    x_new = 0.5*(x_1 + x_2)
    if np.sign(g(x_new)) == sign_x_1:
        x_1 = x_new
    else:
        x_2 = x_new
    error = abs(x_2 - x_1)
#print the result
x_root = 0.5*(x_1+x_2)
print('root using binary search:', x_root)
#more constants
h = 6.626*10**(-34)
c = 3*10**8
k_b = 1.380*10**(-23)
b = h*c/(k_b*x_root)
lambda_sun = 502*10**-9
print("value of Wein's displacement constant:", b)
print('Temperature of the Sun (K):', b/lambda_sun)



