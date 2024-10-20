#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@project: microbiome inference (lotka-volterra model - parameters)
@author: Roman Zapien-Campos - 2023
"""

# Import packages
import pickle as pc
import numpy as np


# Import parameters used in simulations
with open('../data/test_time_rep_noise/simulation_parameters.pickle', 'rb') as f: 
    data_par = pc.load(f)
# Number of microbial types
n_types = data_par['n_types']
# Time simulated
t_simulated = data_par['t_simulated']
# Sampling time points
#t_points = data_par['t_points']
#sampling_times = data_par['sampling_times']
t_points, rep, noise = np.loadtxt('condition.txt',delimiter=',')
t_points = int(t_points)
sampling_times = np.linspace(0, t_simulated, t_points)
# Initial absolute abundance
init_abs_abund = data_par['init_abs_abund']
# Initial relative abundance
init_rel_abund = data_par['init_rel_abund']
# Threshold to stop diverging numerical simulations
upper_threshold = data_par['upper_threshold']