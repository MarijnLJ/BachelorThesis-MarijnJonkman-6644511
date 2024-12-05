#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 11:15:43 2024

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

p_values = np.linspace(0, 1, 22)[1:-1]      # Provides several values (ordered) between 0 and 1, where 0 and 1 are removed
p_values = p_values.tolist()           # Stores the values in a list
expectation_values = [] 
best_q_values = []
a = 2.5       # Minimal size of a side of an outbreak
b = 750     # Maximal size of a side of an outbreak

num_sim = 100 # The amount of simulations that are run

for p in p_values:
    k = 0
    m = 0
    sizes_combined = []  # To store sizes of outbreaks
    not_extinct_count = 0  # Count of non-extinct simulations
    valid_outbreak_sum = 0
    smaller_outbreak_sum = 0
    number_active_sides = 0
    
    for _ in range(num_sim):
        infected, immune = outbreak(p) # Determine the lists of infected and immune nodes
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
      
        # Setting the bin size to 1
        bin_edges = np.arange(np.min(sizes_combined), np.max(sizes_combined) + 2, 1)
        counts, bins = np.histogram(sizes_combined, bins=bin_edges, density=True)
    
        # Plot the histogram with default color
        plt.hist(sizes_combined, bins=bin_edges, density=True, alpha=0.7, color='gray', edgecolor='black')
    
        # Manually adjust the color of each bin
        for i in range(len(bins) - 1):
            bin_center = (bins[i] + bins[i + 1]) / 2
            color = 'lightgray' if bin_center <= (2 * a) else 'dodgerblue'
            plt.bar(bins[i], counts[i], width=bins[i+1] - bins[i], align='edge', color=color, edgecolor='black', alpha=0.7)
    
        plt.title(f'Histogram of Outbreak Sizes for p = {p:.2f}')
        plt.xlabel('Outbreak Size')
        plt.ylabel('Density')
        plt.xlim(left = 0)
        plt.ylim(bottom = 0)
        
        # Define the line that approximates the histograms using lambda
        x_values = np.linspace(0, max(sizes_combined)+1, 100)
        
        # Adding a scaling factor
        scaling_factor = 1 / (np.exp(-lambda_ * a) - np.exp(-lambda_ * b))
        
        line =  scaling_factor * lambda_ * np.exp(-lambda_ * x_values)
        
        plt.plot(x_values, line, color='red', label='Fit Line')
        plt.legend()
        plt.grid(alpha=0.5)
        plt.text(0.85 * plt.xlim()[1], 0.83 * plt.ylim()[1], f'Not extinct: {not_extinct_count}', fontsize=10, ha='center', color='darkred')
        plt.show()
    
