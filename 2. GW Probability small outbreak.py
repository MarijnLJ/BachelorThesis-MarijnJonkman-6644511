#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 10:00:57 2024

@author: marijn
"""

import matplotlib.pyplot as plt
import numpy as np

def root_func(z, p, q):
    """Calculate the root function value."""
    return q - ((1 - p + p * q) ** 3 * np.exp(2 * z * p * (q - 1)))

# Define the ranges for q, p, and zeta values
q_values = np.linspace(0, 1, 1000)[:-1]
p_values = np.linspace(0, 1, 1000)[1:]
zeta_values = [1, 2, 3, 4, 5]
colors = ['blue', 'orange', 'green', 'red', 'purple']

# Loop through each zeta value to calculate the corresponding curves
for i in range(len(zeta_values)):
    z = zeta_values[i]
    
    # Lists to store the p values and the corresponding q values that make root_func zero
    p_for_closest_q = []
    closest_q_values = []

    for p in p_values:
        closest_value = float('inf')  # Start with an infinitely large closest value
        closest_q = None  # Variable to store the q that gives the smallest root_func value
        
        # Loop through all q values to find the one that minimizes the root function
        for q in q_values:
            q_formula = root_func(z, p, q)
            
            # If the current value is closer to zero, update the closest value and the corresponding q
            if abs(q_formula) < closest_value:
                closest_value = abs(q_formula)
                closest_q = q
         
        p_for_closest_q.append(p)
        closest_q_values.append(closest_q)

    # For the found p and q values, we calculate the probability of a small outbreak including the first node
    q_with_index = []
    for p, q in zip(p_for_closest_q, closest_q_values):
        formula = (1 - p + p * q) ** 4 * np.exp(- 2 * z * p * (1 - q))
        q_with_index.append(formula)
    
    # Plotting the results for the current zeta value
    plt.plot(p_for_closest_q, q_with_index, label=f'$\zeta$ = {z}', color = colors[i])

# Formatting of the plot
plt.title('Probability of a small outbreak ')
plt.xlabel('Infection probability ($p$)')
plt.xlim(0, 1)
plt.ylabel(r'Probability of a small outbreak ($\tilde{q}$)')
plt.ylim(0, 1.01)
plt.grid(alpha=0.5)
plt.legend()
plt.show()