# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 23:34:57 2022

PHY407 Lab 5 Q2

@author: Landon Wang
"""

""" The scipy.io.wavfile allows you to read and write .wav files """
from scipy.io.wavfile import read, write
import numpy as np
import matplotlib.pyplot as plt
from numpy import empty


# read the data into two stereo channels
# sample is the sampling rate, data is the data in each channel,
# dimensions [2, nsamples]
sample, data = read('GraviteaTime.wav')
# sample is the sampling frequency, 44100 Hz

# separate into channels
channel_0 = data[:, 0]
channel_1 = data[:, 1]
N_Points = len(channel_0)

# set time variable with number of data points and resize to unit of seconds
t = np.arange(N_Points)/sample

# plot the wave plot
plt.figure(figsize=(8, 6))
plt.subplot(211)
plt.plot(t, channel_0, label = 'Channel_0')
plt.title('Wave Plot for GraviteaTime.wav Seprated for Two Channels')
plt.ylabel('Amplitude')
plt.legend()

plt.subplot(212)
plt.plot(t, channel_1, label = 'Channel_1')
plt.ylabel('Amplitude')
plt.xlabel('Time (s)')
plt.legend()
plt.savefig('Figure 2.1 - Wave Plot.png')

# plot the focused wave plot
plt.figure(figsize=(8, 6))
plt.subplot(211)
plt.plot(t[4*sample:4*sample+1500], channel_0[4*sample:4*sample+1500], label = 'Channel_0')
plt.title('Wave Plot for GraviteaTime.wav Seprated for Two Channels Focused')
plt.ylabel('Amplitude')
plt.legend()

plt.subplot(212)
plt.plot(t[4*sample:4*sample+1500], channel_1[4*sample:4*sample+1500], label = 'Channel_1')
plt.ylabel('Amplitude')
plt.xlabel('Time (s)')
plt.legend()
plt.savefig('Figure 2.2 - Wave Plot Focused.png')


# transformation
z_0=np.fft.fft(channel_0)
z_1=np.fft.fft(channel_1)


# plot the fft plot
plt.figure(figsize=(8, 6))
plt.subplot(211)
plt.plot((t*sample)/(2*np.pi), np.abs(z_0), label = 'Channel_0')
plt.xlim(0, max(t*sample)/(4*np.pi))
plt.title('Fourier Transformated GraviteaTime.wav Seprated for Two Channels')
plt.ylabel('Amplitude')
plt.legend()

plt.subplot(212)
plt.plot((t*sample)/(2*np.pi), np.abs(z_1), label = 'Channel_1')
plt.xlim(0, max(t*sample)/(4*np.pi))
plt.ylabel('Amplitude')
plt.xlabel('Frequency (Hz)')
plt.legend()
plt.savefig('Figure 2.3 - Fourier Transformed wav.png')

# low pass filter
# set the lower limit for low pass filtering Hz
low_lim = 880
# low pass limitor unit change
low_lim = int(low_lim*2*np.pi)

# filterling
z_0[low_lim:-low_lim] = 0
z_1[low_lim:-low_lim] = 0

# plot the filtered fft plot
plt.figure(figsize=(8, 6))
plt.subplot(211)
plt.plot((t*sample)/(2*np.pi), np.abs(z_0), label = 'Channel_0')
plt.xlim(0, 1000)
plt.title('Fourier Transformated GraviteaTime.wav Seprated for Two Channels')
plt.ylabel('Amplitude')
plt.legend()

plt.subplot(212)
plt.plot((t*sample)/(2*np.pi), np.abs(z_1), label = 'Channel_1')
plt.xlim(0, 1000)
plt.ylabel('Amplitude')
plt.xlabel('Frequency (Hz)')
plt.legend()
plt.savefig('Figure 2.4 - Fourier Transformed wav.png')

# transformation back to wave
data_0 = np.fft.ifft(z_0)
data_1 = np.fft.ifft(z_1)

# plot the filtered wave plot
plt.figure(figsize=(8, 6))
plt.subplot(211)
plt.plot(t, data_0, label = 'Channel_0')
plt.title('Filtered Wave Plot for GraviteaTime.wav Seprated for Two Channels')
plt.ylabel('Amplitude')
plt.legend()

plt.subplot(212)
plt.plot(t, data_1, label = 'Channel_1')
plt.ylabel('Amplitude')
plt.xlabel('Time (s)')
plt.legend()
plt.savefig('Figure 2.5 - Filtered Wave Plot.png')


# plot the filtered focused wave plot
plt.figure(figsize=(8, 6))
plt.subplot(211)
plt.plot(t[4*sample:4*sample+1500], data_0[4*sample:4*sample+1500], label = 'Channel_0')
plt.title('Filtered Wave Plot for GraviteaTime.wav Seprated for Two Channels Focused')
plt.ylabel('Amplitude')
plt.legend()

plt.subplot(212)
plt.plot(t[4*sample:4*sample+1500], data_1[4*sample:4*sample+1500], label = 'Channel_1')
plt.ylabel('Amplitude')
plt.xlabel('Time (s)')
plt.legend()
plt.savefig('Figure 2.6 - Filtered Wave Plot Focused.png')

# this creates an empty array data_out with the same shape as "data"
# (2 x N_Points) and the same type as "data" (int16)
data_out = empty(data.shape, dtype = data.dtype)
# fill data_out
data_out[:, 0] = data_0
data_out[:, 1] = data_1
write('GraviteaTime_lpf.wav', sample, data_out)



