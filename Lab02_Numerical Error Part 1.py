'''Written by Yinshi Liu'''
import numpy as np
import scipy as sp

#create an array from loaded data
a = np.random.normal(0., 1., 2000)
b = np.random.normal(1.e7, 1., 2000)

#create an array from loaded data
data = np.loadtxt('cdata.txt')

#Calculate average
def average(data):
    return np.sum(data)/len(data)

#Calculate std using equation 1 
def equation_1 (data, avg):
    '''Same code as Q1b, modularized to individual functions'''
    sum_1 = 0
    for i in data:
        sum_1 += (i - avg)**2
    return np.sqrt(1/(len(data)-1) * sum_1)

#calculate std using equation 2
def equation_2 (data):
    '''Same code as Q1b, modularized to individual functions'''
    sum_2 = 0
    sum_1 = 0
    for i in data:
        sum_1 += i
        sum_2 += i**2
    var_2 = 1/(len(data)-1)*(sum_2 - sum_1**2/len(data))
    if var_2 < 0:
        print('negative variance encountered! Var_2 =', var_2)
        return 0
    else:
        return np.sqrt(var_2)

#function for std
def standard_deviation (data):
    '''Calls the np.std function'''
    return np.std(data, ddof = 1)

print('result from equation 1 =', equation_1 (data, average(data)))
print('result from euqation 2 =', equation_2 (data))
print('true value =', standard_deviation(data))
print('standard error for eqn. 1 =',
      (equation_1(data,average(data))-standard_deviation(data))/standard_deviation(data))
print('standard error for eqn. 2 =',
      (equation_2(data)-standard_deviation(data))/standard_deviation(data))
print ('')

sigma_1a = equation_1 (a, average(a))
std_a = standard_deviation(a)
r_error_1a = (sigma_1a - std_a)/std_a

sigma_1b = equation_1 (b, average(b))
std_b = standard_deviation(b)
r_error_1b = (sigma_1b - std_b)/std_b

print('sigma_1a =', sigma_1a, 'np.std for a =', std_a)
print('magnitude of relative error with eqn 1 for a =', abs(r_error_1a))

sigma_2a = equation_2 (a)
r_error_2a = (sigma_2a - std_a)/std_a

sigma_2b = equation_2 (b)
r_error_2b = (sigma_2b - std_b)/std_b

print('sigma_2a =', sigma_2a, 'np.std =', std_a)
print('magnitude of relative error with eqn 2 for a', abs(r_error_2a))
print('')

print('sigma_1b =', sigma_1b, 'np.std for a =', std_b)
print('magnitude of relative error with eqn 1 for b =', abs(r_error_1b))
print('sigma_2b =', sigma_2b, 'np.std for a =', std_b)
print('magnitude of relative error with eqn 1 for b =', abs(r_error_2b))
print('')
