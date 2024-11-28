#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 11:13:30 2024

@author: marijn
"""

import matplotlib.pyplot as plt
import numpy as np

q_values = np.linspace(0, 1, 100) 

def gen_function(z, p, q):
    """ Calculate the probability generating function (pgf) """
    return (1-p+p*q)**3*np.exp(2*z*p*(q-1))

# Represents the line y = q 
simple_function = [q for q in q_values]

# Plot the pgf for R_0>1
plt.plot(q_values, gen_function(17/2 ,0.1, q_values))
plt.plot(q_values, simple_function, linestyle = '--')
plt.xlabel('q')
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.title('Probability generating function for $R_0>1$')
plt.grid(alpha = 0.5)
plt.show()

# Plot the pgf for R_0<1
plt.plot(q_values, gen_function(2, 0.1, q_values))
plt.plot(q_values, simple_function, linestyle = '--')
plt.xlabel('q')
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.title('Probability generating function for $R_0<1$')
plt.grid(alpha = 0.5)
plt.show()

p_values = np.linspace(0, 1, 1000)[1:]

# Calculate critical equation values, setting values below zero to zero
crit_eq = [max(1/(2*p) - 3/2, 0) for p in p_values]

# Plot the critical equation for values of p (infection probability)
plt.text(0.45, 25, '$R_0>1$', fontsize = 18)
plt.text(0.01, 2, '$R_0<1$', fontsize = 10)
plt.plot(p_values, crit_eq, linewidth = 2, color = 'blue')
plt.xlabel('Infection probability (p)')
plt.ylabel('Expected number of added edges ($\zeta$)')
plt.xlim(0, 1)
plt.ylim(0, 50)
plt.title('Expected number of added edges')
plt.grid(alpha = 0.5)
plt.show()
