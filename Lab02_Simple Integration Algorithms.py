'''Written by Yinshi Liu'''
import numpy as np
from time import time
#Defining functions
def f(x):
    return 4/(1+x**2)

def trapezoidal (N, a, b):
    '''Trapezoidal method adapted from the code fragment in textbook'''
    h = (b-a)/N

    s = 0.5*f(a) + 0.5*f(b)
    for i in range(1, N):
        s += f(a+i*h)
    return h*s

def simpson (N, a, b):
    '''Perform the Simpson's method of intergration.
       Method adopted from formula on textbook'''
    h = (b-a)/N
    i = 1
    even = 0
    odd = 0
    while i in range(1, N):
        #even number
        if i % 2 == 0:
            even += f(a+i*h)
        #odd number
        else:
            odd += f(a+i*h)
        i += 1
    return (f(a) + f(b) + 4*odd + 2*even) * h/3

#Part b
print('trapezoidal method =', trapezoidal(4, 0, 1))
print("Simpson's Method =", simpson(4, 0, 1))
print('true value =', np.pi)
print('')
#Part c
#Define initial conditions
n = 1
trapeziodal_error = trapezoidal (2**n, 0, 1) - np.pi
#timing trapezoidal method

#initialize loop, start timing
start_trapeziodal = time()

#a simple while loop, since this method underestimates the value, I added
#an abs() function
while abs(trapeziodal_error) > 1*10**-9:
    n += 1
    trapeziodal_error = trapezoidal (2**n, 0, 1) - np.pi

#end loop, stop timing
end_trapezoidal = time()

#print the results
print('results for trapezoidal method: ')
print('number of slices = 2^',n,'slices')
print('time taken =',end_trapezoidal - start_trapeziodal,'s')
print('result =', trapezoidal(2**n, 0, 1))

#timing simpson's method
#initialize loop
start_simpson = time()
n_1 = 0
simpson_error = simpson (2**n_1, 0, 1) - np.pi
#while loop for timing
while abs(simpson_error) > 1*10**-9:
    n_1 += 1
    simpson_error = simpson (2**n_1, 0, 1) - np.pi
#stop timing
end_simpson = time()

#print the results
print('')
print("results for Simpson's method: ")
print('number of slices = 2^',n_1,'slices')
print('time taken =',end_simpson - start_simpson,'s')
print('result =', simpson(2**n, 0, 1))
print('')

#part d and e
#Simpson's
I_1 = simpson(16, 0, 1)
I_2 = simpson(32, 0, 1)
error = 1/15*abs(I_2 - I_1)

#trapezoidal
I_1t = trapezoidal(16, 0, 1)
I_2t = trapezoidal(32, 0, 1)
error_t = 1/3*abs(I_2t - I_1t)

#print the results
print ('estimated error for trap. =', error_t)
print ('actual error for trap. =', np.pi-I_2t)
print ('estimated error for Simpson=', error)
print ('actual error for Simpson=', np.pi - I_2)


