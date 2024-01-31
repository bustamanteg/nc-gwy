#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 22:46:02 2022

@author: lfairgrievepark12
"""

import numpy as np
import pandas as pd
import glob
import os
import sys
import matplotlib.pyplot as plt
import scipy.optimize as sci

def parabola(x,a,b,c):
    y = a*x**2+b*x+c
    return y

data = np.genfromtxt('Voltagesweep_-5-5V_QDsite1031-VP.vpdata.csv', delimiter=',', skip_header=1)

bias = data[:,1]
freq = data[:,2]/-0.02
diss = data[:,3]
#diss = (diss-diss.min())/diss.min()

popt, pcovy = sci.curve_fit(parabola,bias,freq)
freq_sub = freq - parabola(bias, *popt)

plt.figure(1)
plt.plot(bias)
#plt.xlim([-5,-2])
plt.xlabel('Bias (V)')
plt.ylabel('Freq shift (Hz)')

plt.figure(2)
plt.plot(bias,freq_sub)
#plt.xlim([-5,-2])
plt.xlabel('Bias (V)')
plt.ylabel('Freq shift_sub (Hz)')

plt.figure(3)
plt.plot(bias,diss)
#plt.xlim([-5,-2])
plt.xlabel('Bias (V)')
plt.ylabel('Diss (Hz)')