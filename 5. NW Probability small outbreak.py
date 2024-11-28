#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 11:02:03 2024

@author: marijn
"""

import random                                                     
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import math
import scipy.odr as odr

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

p_values = np.linspace(0, 1, 82)[1:-1]      # Provides several values (ordered) between 0 and 1, where 0 and 1 are removed
p_values = p_values.tolist()           # Stores the values in a list 
a = 2.5       # Minimal size of a side of an outbreak
b = 750     # Maximal size of a side of an outbreak

num_sim = 100
zeta = [0.1, 0.5, 1, 2]
colors = ['blue', 'orange', 'green', 'red']

# Loop through each zeta value to calculate the corresponding curves
for i in range(len(zeta)):
    z = zeta[i]
    
    p_for_best_q = []
    best_q_values = []
    
    for p in p_values:
        k = 0
        m = 0
        sizes_combined = []  # Stores sizes of outbreaks
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
        
        valid_outbreaks = 2 * num_sim - m - k     # The number of sides that are interesting for our analysis
            
        outbreak_counts = Counter(sizes_combined)  # Count the frequency of each outbreak size
        total_outbreaks = len(sizes_combined)  # Total number of outbreaks simulated
        
        # Calculate the outbreak probabilities based on the outbreak counts
        outbreak_probabilities = {size: count / total_outbreaks for size, count in outbreak_counts.items()}
        
        prob_Y_list = []
        outbreak_sizes = sorted(outbreak_probabilities.keys()) # Sort the outbreak sizes
        threshold = 1e-10 # Define a small threshold value to filter out very small probabilities

        for y in range(51):
            prob_Y = 0
            for size in outbreak_sizes:
                probability = outbreak_probabilities.get(size, 0) # Get the probability of this outbreak size (default to 0 if not found)
                
                # Add the contribution of the current outbreak size to prob_Y
                prob_Y += (1 / math.factorial(y)) * probability * np.exp(-2 * z * p * size) * ((2 * z * p * size) ** y)
            
            # If prob_Y is above the threshold, store the (y, prob_Y) pair
            if prob_Y > threshold:
                prob_Y_list.append((y, prob_Y))
        
        q_values = np.linspace(0, 1, 10000)[:-1]
        min_diff = float('inf')
        best_q = None
        
        for q in q_values:
            # Calculate the q formula by summing over all y values
            q_formula = sum((q ** y) * prob_Y for y, prob_Y in prob_Y_list) - q
            diff = abs(q_formula)
            
            # If the current difference is smaller than the previous minimum, update min_diff and best_q
            if diff < min_diff:
                min_diff = diff
                best_q = q
        
        # Store the current value of p and the q values for which the formula is closest to zero
        p_for_best_q.append(p)
        best_q_values.append(best_q)
    
    # Plot the values we found for the current zeta
    plt.plot(p_for_best_q, best_q_values, 'o', markersize = 4, label=f'$\zeta$ = {z}', color = colors[i])
    
    # To fit a line through the points we need to ignore some of the values. 
    # For small p we note that q is almost 1, so we want to start the fit when q starts to decrease
    for j in range(len(best_q_values)):
        if best_q_values[j + 1] < best_q_values[j]:
            mark = j
            break
    
    # Store the values we want for the fit in two lists
    p_for_best_q2 = p_for_best_q[mark:]
    best_q_values2 = best_q_values[mark:]
    
    # Starting values for the fit
    # if z == 0.1:
    #     B0_start= 2e-1
    #     B1_start= -1200
    #     B2_start = 5200
    # else:
    B0_start= 1e-4
    B1_start= -10
    B2_start = 5
    
    def f(B, x):
        return B[0] + B[2]*np.exp(B[1]*x)

    odr_model = odr.Model(f) # Define the model-object to use in odr
    odr_data  = odr.RealData(p_for_best_q2, best_q_values2) # Define a RealData object
    odr_obj = odr.ODR(odr_data, odr_model, beta0=[B0_start, B1_start, B2_start]) # Make a odr object with the data, model and starting values
    odr_res = odr_obj.run() # Execut the fit

    par_best = odr_res.beta # Take the best estimator for the parameter from the data
    print(par_best)
    # Plot the fitted line
    xplot = np.arange(0, 1, 0.001)
    plt.plot(xplot, f(par_best, xplot), color = colors[i], linestyle = 'dashed')


# Formatting of the plot    
plt.xlabel('Infection Probability ($p$)')
plt.ylabel('Probability of a small outbreak ($q$)')
plt.title('Probability of a small outbreak')
plt.xlim(0, 1)
plt.ylim(0, 1.01)
plt.grid(alpha = 0.5)
plt.legend()
plt.show()