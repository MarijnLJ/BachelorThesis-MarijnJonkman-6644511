#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 11:02:03 2024

@author: marijn
"""

import random                                                     
import matplotlib.pyplot as plt
import numpy as np

# Set a seed for random number generation to ensure reproducibility
random.seed(3)          
                                                             
def outbreak(p):
    """
    This function determines the lists infected and immune for several generations. 
    input:
        p: infection probability (float)
    output:
        infected: individuals that are infected (list)
        immune: individuals that are immune (list)

    """
    infected = [0]                      # Initially, only person 0 is infected
    immune = []                         # No one is immune at the start
    generation = 0
    while infected and generation < 400:
        infected_new = []
        for x in infected:              # Iterates through all infected individuals
            for y in [x-2, x-1, x+1, x+2]:          # Check the neighbours of x
                if y not in infected and y not in immune and y not in infected_new:     # Check if y is susceptible
                    if random.random() <= p:        # With probability p, y gets infected
                        infected_new.append(y)
        immune.extend(infected)                     # Infected individuals become immune in the next generation
        infected = infected_new                     # Newly infected individuals become the infected group for the next generation 
        generation += 1                           
    return infected, immune

p_values = np.linspace(0, 1, 102)[1:-1]      # Provides several values (ordered) between 0 and 1, where 0 and 1 are removed
p_values = p_values.tolist()           # Stores the values in a list
expectation_values = [] 
best_q_values = []
a = 2.5       # Minimal size of a side of an outbreak
b = 750     # Maximal size of a side of an outbreak

num_sim = 100

for p in p_values:
    k = 0
    m = 0
    sizes_combined = []  # To store sizes of outbreaks
    not_extinct_count = 0  # Count of non-extinct simulations
    valid_outbreak_sum = 0
    smaller_outbreak_sum = 0
    number_active_sides = 0
    
    for _ in range(num_sim):
        infected, immune = outbreak(p)
        total_outbreak = infected + immune
        
        # Separate the outbreak into negative (left) and positive (right) individuals
        left_side = [individual for individual in total_outbreak if individual < 0]
        right_side = [individual for individual in total_outbreak if individual >= 0]
        
        outbreak_size_left = len(left_side) # Size of the left side
        outbreak_size_right = len(right_side) # Size of the right side
        
        # Calculate the sum of all sides that are smaller than or equal to a
        if outbreak_size_left <= a:
            k += 1
        if outbreak_size_right <= a:
            k += 1
        
        # Calculate the sum of all sides that are bigger than b
        if outbreak_size_left > b:
            m += 1
        if outbreak_size_right > b:
            m += 1
        
        # Calculate the sum of all sides of an oubreak between a and b
        if outbreak_size_left > a and outbreak_size_left <= b:
            valid_outbreak_sum += outbreak_size_left
        if outbreak_size_right > a and outbreak_size_right <= b:
            valid_outbreak_sum += outbreak_size_right
        
        # Calculate the sum of all sides of an oubreak smaller than b
        if outbreak_size_left <= b:
            smaller_outbreak_sum += outbreak_size_left
        if outbreak_size_right <= b:
            smaller_outbreak_sum += outbreak_size_right
        
        outbreak_size = outbreak_size_left + outbreak_size_right
        
        # If the outbreak is not extinct we up the count and do not add the outbreak size to sizes_combined
        if len(infected) != 0:
            not_extinct_count += 1
            continue
        
        sizes_combined.append(outbreak_size)
    
    valid_outbreaks = 2 * num_sim - m - k     # number of sides that are interesting for our analysis
    
    # If there are valid outbreaks, we calculate the estimator: lambda. 
    # If there are no valid outbreaks, the expected outbreak size is calculated without lambda
    if valid_outbreaks != 0:
        lambda_ = valid_outbreaks / (valid_outbreak_sum + m * b - a * (2 * num_sim - k))
    
    # For this lambda, the expected outbreak size is calculated
        expect = (smaller_outbreak_sum + m * (b + 1 / lambda_)) / (2 * num_sim)
        expectation_values.append((p, expect))
    else:
        expect = smaller_outbreak_sum / (2 * num_sim)
        expectation_values.append((p, expect))

# Find the values of p and the expectation value to use in the plots
plot_p_vals = [item[0] for item in expectation_values]
plot_expect_vals = [item[1] for item in expectation_values]

# Calculate zeta for these values
zeta_values = []
for e, p in zip(plot_expect_vals, plot_p_vals):
    if e != 0:
        zeta = 1 / (2 * e * p) 
        zeta_values.append(zeta)
    else:
        zeta = 0
        zeta_values.append(zeta)

# Plot the expected outbreak size against p
plt.plot(plot_p_vals, plot_expect_vals, color='b')
plt.xlabel("Infection Probability ($p$)", fontsize=12)
plt.ylabel("Expected Outbreak Size", fontsize=12)
plt.xlim(0, 1)
plt.ylim(0, 100000)
plt.title("Expected outbreak size", fontsize=14)
plt.grid(alpha = 0.5)
plt.show()

# To see the graph better for smaller values of p we filter this out
filtered_p_vals = [p for p in plot_p_vals if p < 0.95]
filtered_expect_vals = [expect for p, expect in zip(plot_p_vals, plot_expect_vals) if p < 0.95]

plt.plot(filtered_p_vals, filtered_expect_vals, color='b')
plt.xlabel("Infection Probability ($p$)", fontsize=12)
plt.ylabel("Expected Outbreak Size", fontsize=12)
plt.xlim(0, 1)
plt.ylim(0, 1000)
plt.title("Expected outbreak size", fontsize=14)
plt.grid(alpha = 0.5)
plt.show()

# Plot zeta against p
plt.plot(plot_p_vals, zeta_values, color = 'b', linewidth = 2)                     
plt.text(0.45, 25, '$R_0>1$', fontsize = 18)
plt.text(0.015, 2, '$R_0<1$', fontsize = 12)
plt.xlabel('Infection Probability (p)')
plt.ylabel('Expected number of shortcuts ($\zeta$)')
plt.xlim(0, 1)
plt.ylim(0, 50)
plt.title('Expected number of shortcuts')
plt.grid(alpha = 0.5)
plt.show()
