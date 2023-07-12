'''Writte by Yinshi Liu'''
import numpy as np
import scipy as sp
from numpy import array, empty
from time import time
import matplotlib.pyplot as plt
from SolveLinear import GaussElim
#part a
#Check matrix
A = array([[2, 1, 4, 1],
           [3, 4, -1, -1],
           [1, -4, 1, 5],
           [2, -2, 1, 3]], float)
v = array([-4, 3, 9, 7], float)
#define functions
def partial_pivot (A, v):
    '''Adapted from Newman'''
    N = len(v)
    for m in range(N):
        diag = abs(A[m,m])
        index = m
    #Check if pivot is small (arbitary number)
        if diag < 0.1:
            #find largest value
            for i in range(m+1, N):
                if abs(A[i,m]) > diag:
                    diag = abs(A[i,m])
                    #Swap rows
                    if i != m:
                        A[i], A[m] = (A[m].copy(),A[i].copy())
                        v[m], v[i] = (v[i], v[m])
                        break
    #Continue elimination               
        div = A[m,m]
        A[m,:] /= div
        v[m] /= div

        for i in range(m+1, N):
            mult = A[i, m]
            A[i,:] -= mult*A[m,:]
            v[i] -= mult*v[m]
    x = empty(N, float)
    for m in range(N-1, -1, -1):
        x[m] = v[m]
        for i in range(m+1, N):
            x[m] -= A[m, i]*x[i]
    return x

def gaussian_elimination (A, v):
    '''From Sample Code'''
    return GaussElim(A, v)
    

def LU_decomposition (A, v):
    '''LU decomposes matrix then solves matrix'''
    N = len(v)
    #Define Upper and Lower matrices 
    #Use elimination method to find U, copy A to eliminate
    U = A
    #Construct L using algorithm outlined in eqn 6.32
    L = np.identity(N)
    for m in range(N):
        div = U[m,m]
        for j in range(m+1, N):
            L[j,m] = U[j,m]/div
            U[j] = U[j] - L[j,m]*U[m]
    #calculate y values
    y = empty(N, float)
    #Initial condition for y[0], the first value of y
    y[0] = v[0]/L[0,0]
    #forward method
    for k in range(1,N):
        y[k] = (v[k] - np.dot(L[k, 0:k], y[0:k]))/L[k,k]
    #calculate x values
    x = empty(N, float)
    #Initial condition for x[N-1], the last value of x
    x[N-1] = y[N-1]/U[N-1,N-1]
    #backward method
    for k in range(N-2, -1, -1):
        x[k] = (y[k] - np.dot(U[k, k+1:N], x[k+1:N]))/U[k,k]
    return x
#Code outputs
print('Solution using Gaussian elimination =', gaussian_elimination (A, v))
print('Solution using Partial Pivot =', partial_pivot (A, v))
print('Solution using Numpy =', np.linalg.solve(A, v))
#part b
#Initialize arrays
size = array([10, 20, 30])
time_pp = []
time_ge = []
time_LU = []
err_pp = []
err_ge = []
err_LU = []

#Size of matrix
size = array([5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
              150, 200, 250, 300, 350, 400, 450, 500, 600])
#initialize loop
for i in size:
    #Create random positive integer matrix w/out zeros
    matrix = np.random.randint(1, 10, size = (i, i))
    vector = np.random.randint(1, 10, size = i)
    #convert into float
    matrix = matrix.astype(np.float)
    vector = vector.astype(np.float)

    #time and error for partial pivot method
    start_pp = time()
    x_pp = partial_pivot(matrix, vector)
    end_pp = time()
    time_pp.append(end_pp - start_pp)
    #error
    vsol_pp = np.dot(matrix, x_pp)
    err_pp.append(np.mean(abs(vector - vsol_pp)))
    #time error gaussian
    start_ge = time()
    x_ge = gaussian_elimination(matrix, vector)
    end_ge = time()
    time_ge.append(end_ge - start_ge)
    
    vsol_ge = np.dot(matrix, x_ge)
    err_ge.append(np.mean(abs(vector - vsol_ge)))
    #time error LU
    start_LU = time()
    x_LU = LU_decomposition(matrix, vector)
    end_LU = time()
    time_LU.append(end_LU - start_LU)
    vsol_LU = np.dot(matrix, x_LU)
    err_LU.append(np.mean(abs(vector - vsol_LU)))

time_pp = array(time_pp)
time_ge = array(time_ge)
time_Lu = array(time_LU)
err_pp = array(err_pp)
err_ge = array(err_ge)
err_LU = array(err_LU)
#plot the graph
plt.figure()
plt.scatter(size, time_ge,
            label = 'Time using Gaussian elimination', s = 5)
plt.scatter(size, time_pp,
            label = 'Time using Partial pivot', s = 5)
plt.scatter(size, time_LU,
            label = 'Time using LU decomposition', s = 5)
#plt.yscale('log')
#plt.xscale('log')
plt.xlabel('Matrix Size (N)')
plt.ylabel( 'Time (s)')
plt.legend()
plt.title('Computing Time for Diffenent Matrix Solving Algorithms')
plt.savefig('Figure 1.1 Computing Time for Solving Matrices.png')

plt.figure()
plt.scatter(size, err_ge,
            label = 'Error using Gaussian elimination', s = 5)
plt.scatter(size, err_pp,
            label = 'Error using Partial pivot', s = 5)
plt.scatter(size, err_LU,
            label = 'Error using LU decomposition', s = 5)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Matrix Size (N)')
plt.ylabel( 'Error Size')
plt.legend()
plt.title('Error size for Diffenent Matrix Solving Algorithms')
plt.savefig('Figure 1.2 Error Size for Solving Matrices.png')
