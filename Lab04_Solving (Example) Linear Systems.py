'''Written by Yinshi Liu, code adapted from Newman'''
import numpy as np
import scipy as sp
from numpy import array, empty

A = array([[0, 7,-1, 3, 1],
           [0, 3, 4, 1, 7],
           [6, 2, 0, 2,-1],
           [2, 1, 2, 0, 2],
           [3, 4, 1,-2, 1]], float)
v = array([5, 7, 2, 3, 4], float)
N = len(v)
print(np.linalg.solve(A, v))
#Partial Pivot
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
print(x)
