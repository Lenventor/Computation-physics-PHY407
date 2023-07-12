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

#plot for part e
#Define constants
x_0 = 1
x_c = 8.66*10**7

#initialize the loop
x_values = [x_0, x_c, 10*x_c]
T_values = []
#while loop to iterate through x values
for x in x_values:
    N = 16
    a = 0.0
    b = x
    T = 0
    x, w = gaussxw(N)
    xp = 0.5*(b-a)*x + 0.5*(b+a)
    wp = 0.5*(b-a)*w
    for i in range(N):
        T += wp[i]*t(xp[i])
    T_values.append(T)
    
T_values = np.array(T_values)

print(T_values)

