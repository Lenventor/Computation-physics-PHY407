# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 02:16:11 2022

 PHY407 Lab 5 Q3

@author: Landon Wang
"""

import numpy as np
import matplotlib.pyplot as plt

# data imputs
SLP = np.loadtxt('SLP.txt')
Longitude = np.loadtxt('lon.txt')
Times = np.loadtxt('times.txt')

# plot SLP
plt.figure(figsize=(8, 6))
plt.contourf(Longitude, Times, SLP)
plt.xlabel('longitude(degrees)')
plt.ylabel('days since Jan. 1 2015')
plt.title('SLP anomaly (hPa)')
plt.colorbar()
plt.savefig('Figure 3.1 - Oraginal SLP.png')
plt.show()


# Filter
# create data frame
SLP_f3 = np.zeros((120,144))
SLP_f5 = np.zeros((120,144))
# loop every days in the data
for d in Times:
    # extract data
    SLP_d = SLP[int(d)]
    # fourier transformation
    z_d = np.fft.fft(SLP_d)
    # filter only m = 3 and 5
    z_d_f3 = np.zeros(len(z_d), dtype=np.complex_)
    z_d_f3[3] = z_d[3]
    z_d_f3[-3] = z_d[-3]
    
    z_d_f5 = np.zeros(len(z_d), dtype=np.complex_)
    z_d_f5[5] = z_d[5]
    z_d_f5[-5] = z_d[-5]
    # reverse fourier transformation
    SLP_d_f3=np.fft.ifft(z_d_f3)
    SLP_d_f5=np.fft.ifft(z_d_f5)
    # store data
    SLP_f3[int(d)] = SLP_d_f3
    SLP_f5[int(d)] = SLP_d_f5


# plot SLP for m = 3
plt.figure(figsize=(8, 6))
plt.contourf(Longitude, Times, SLP_f3)
plt.xlabel('longitude(degrees)')
plt.ylabel('days since Jan. 1 2015')
plt.title('Filtered SLP anomaly to m = 3 (hPa)')
plt.colorbar()
plt.savefig('Figure 3.2 - Filtered SLP M = 3.png')
plt.show()

# plot SLP for m = 5
plt.figure(figsize=(8, 6))
plt.contourf(Longitude, Times, SLP_f5)
plt.xlabel('longitude(degrees)')
plt.ylabel('days since Jan. 1 2015')
plt.title('Filtered SLP anomaly to m = 5 (hPa)')
plt.colorbar()
plt.savefig('Figure 3.3 - Filtered SLP M = 5.png')
plt.show()





