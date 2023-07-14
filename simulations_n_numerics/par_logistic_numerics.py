#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@project: microbiome inference (logistic model - parameters)
@author: Roman Zapien-Campos - 2023
"""

# Import packages
import pickle as pc

### PLEASE UNCOMMENT ONE OF THE NEXT TWO BLOCKS WHILE COMMENTING THE OTHER ONE ###

## UNCOMMENT TO RUN NOTEBOOK 02_logistic.ipynb ##

# Import parameters used in simulations
with open('../data/logistic/simulation_parameters.pickle', 'rb') as f:
    data_par = pc.load(f)
# Number of microbial types
n_types = data_par['n_types']
# Time simulated
t_simulated = data_par['t_simulated']
# Sampling time points
t_points = data_par['t_points']
sampling_times = data_par['sampling_times']
# Initial absolute abundance
init_abs_abund = data_par['init_abs_abund']
# Initial relative abundance
init_rel_abund = data_par['init_rel_abund']
# Threshold to stop diverging numerical simulations
upper_threshold = data_par['upper_threshold']


## UNCOMMENT TO RUN NOTEBOOK 05_omm12.ipynb ##

# # Import parameters from empirical data
# with open('../data/omm12/experimental_abs_abund.pickle', 'rb') as f:
#     data_par = pc.load(f)
# # Number of microbial types
# n_types = data_par['n_types']
# # Sampling time points
# t_simulated = data_par['sampling_times'][-1]
# t_points = data_par['t_points']
# sampling_times = data_par['sampling_times']
# # Initial absolute abundance
# init_abs_abund = data_par['init_abs']#[:n_types]
# # Initial relative abundance
# with open('../data/omm12/experimental_rel_abund.pickle', 'rb') as f:
#     data_par = pc.load(f)
# init_rel_abund = data_par['init_comp']
# # Threshold to stop diverging numerical simulations
# upper_threshold = 3E7