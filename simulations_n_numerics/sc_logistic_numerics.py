#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@project: microbiome inference (logistic model - source code)
@author: Roman Zapien-Campos - 2023
"""

# Import packages
import numpy as np
from scipy import integrate
from par_logistic_numerics import *


# Indexes of elements in a matrix of size n_types
diag = np.eye(n_types)
off_diag = 1 - diag


## Moment Closure - 3rd moments (WORKING)
def eq_moments_abs_abund_closure_from_3rd(t, y, gR, dR, mR, N):
    # Unpack different moments        
    m_k = y[:n_types]
    cm_kj = y[n_types:].reshape(n_types, n_types)
    m_kk = np.diag(cm_kj)
    
    # First moment for each type
    D_m_k = gR * (m_k - cm_kj.sum(1) / N)
    D_m_k += mR * (1. - m_k.sum() / N)
    D_m_k -= dR * m_k / N
    
    # Second moment for each type
    D_cm_kj = diag * gR * (m_k - cm_kj.sum(1) / N + 2. * m_kk * (1. - m_k.sum() / N))
    D_cm_kj += diag * mR * (1. - m_k.sum() / N + 2. * (m_k - cm_kj.sum(1) / N))
    D_cm_kj += diag * dR * (m_k - 2. * m_kk) / N
    
    # Co-moment for types k and l
    for k in range(n_types):
        for l in range(n_types):
            if k != l:
            
                D_cm_kj[k,l] += (gR[k] + gR[l]) * cm_kj[k,l] * (1. - m_k.sum() / N)
                D_cm_kj[k,l] += mR[k] * (m_k[l] - cm_kj[l,:].sum() / N) + mR[l] * (m_k[k] - cm_kj[k,:].sum() / N)
                D_cm_kj[k,l] -= (dR[k] + dR[l]) * cm_kj[k,l] / N

    return np.hstack((D_m_k, D_cm_kj.reshape(1, n_types**2)[0]))


## Moment Closure - 3rd moments (WORKING)
def eq_moments_rel_abund_closure_from_3rd(t, y, gR, dR, mR, N):

    # Unpack different moments
    m_k = y[:n_types]
    cm_kj = y[n_types:-1].reshape(n_types, n_types)
    m_kk = np.diag(cm_kj)
    m_Sigma = y[-1]
   
    # First moment for each type
    D_m_k = gR * (m_k * m_Sigma - cm_kj.sum(1) * m_Sigma**2 / N)
    D_m_k += mR * (1. - m_k.sum() * m_Sigma / N)
    D_m_k -= dR * m_k * m_Sigma / N
    
    # First moment of scaling factor
    D_m_Sigma = sum(D_m_k)
    
    # First moment for each type (continuation)
    D_m_k -= m_k * D_m_Sigma
    D_m_k /= m_Sigma
    
    # Second moment for each type
    D_cm_kj = diag * gR * (m_k * m_Sigma - cm_kj.sum(1) * m_Sigma**2 / N + 2. * m_kk * m_Sigma**2 * (1. - m_k.sum() * m_Sigma / N))
    D_cm_kj += diag * mR * (1. - m_k.sum() * m_Sigma / N + 2. * (m_k * m_Sigma - cm_kj.sum(1) * m_Sigma**2 / N))
    D_cm_kj += diag * dR * (m_k * m_Sigma - 2. * m_kk * m_Sigma**2) / N
    
    # Co-moment for types k and l
    for k in range(n_types):
        for l in range(n_types):
            if k != l:
            
                D_cm_kj[k,l] += (gR[k] + gR[l]) * cm_kj[k,l] * m_Sigma**2 * (1. - m_k.sum() * m_Sigma / N)
                D_cm_kj[k,l] += mR[k] * (m_k[l] * m_Sigma - cm_kj[l,:].sum() * m_Sigma**2 / N) + mR[l] * (m_k[k] * m_Sigma - cm_kj[k,:].sum() * m_Sigma**2 / N)
                D_cm_kj[k,l] -= (dR[k] + dR[l]) * cm_kj[k,l] * m_Sigma**2 / N
    
    # Second moment for each type (continuation)
    D_cm_kj -= diag * 2. * m_kk * m_Sigma * D_m_Sigma
    
    # Co-moment for types k and j (continuation)
    D_cm_kj -= off_diag * 2. * cm_kj * m_Sigma * D_m_Sigma
    
    D_cm_kj /= m_Sigma**2

    return np.hstack((D_m_k, D_cm_kj.reshape(1, n_types**2)[0], D_m_Sigma))


## Moment Closure - 2nd & 3rd moments (WORKING)
def eq_moments_abs_abund_closure_from_2nd(t, y, gR, dR, mR, N):

    # Unpack different moments
    m_k = y[:n_types]
    m_k_prod = np.outer(m_k, m_k)
    
    # First moment for each type
    D_m_k = gR * (m_k - m_k_prod.sum(1) / N)
    D_m_k += mR * (1. - m_k.sum() / N)
    D_m_k -= dR * m_k / N
    
    # Second moment for each type
    D_cm_kj = diag * gR * (m_k - m_k_prod.sum(1) / N + 2. * m_k**2 * (1. - m_k.sum() / N))
    D_cm_kj += diag * mR * (1. - m_k.sum() / N + 2. * (m_k - m_k_prod.sum(1) / N))
    D_cm_kj += diag * dR * (m_k - 2. * m_k**2) / N
    
    # Co-moment for types k and l
    for k in range(n_types):
        for l in range(n_types):
            if k != l:
            
                D_cm_kj[k,l] += (gR[k] + gR[l]) * m_k_prod[k,l] * (1. - m_k.sum() / N)
                D_cm_kj[k,l] += mR[k] * (m_k[l] - m_k_prod[l,:].sum() / N) + mR[l] * (m_k[k] - m_k_prod[k,:].sum() / N)
                D_cm_kj[k,l] -= (dR[k] + dR[l]) * m_k_prod[k,l] / N

    return np.hstack((D_m_k, D_cm_kj.reshape(1, n_types**2)[0]))


## Moment Closure - 2nd & 3rd moments (WORKING)
def eq_moments_rel_abund_closure_from_2nd(t, y, gR, dR, mR, N):
   # Unpack different moments         
    m_k = y[:n_types]
    m_k_prod = np.outer(m_k, m_k)
    m_Sigma = y[-1]
   
    # First moment for each type
    D_m_k = gR * (m_k * m_Sigma - m_k_prod.sum(1) * m_Sigma**2 / N)
    D_m_k += mR * (1. - m_k.sum() * m_Sigma / N)
    D_m_k -= dR * m_k * m_Sigma / N
    
    # First moment of scaling factor
    D_m_Sigma = sum(D_m_k)
    
    # First moment for each type (continuation)
    D_m_k -= m_k * D_m_Sigma
    
    D_m_k /= m_Sigma
    
    # Second moment for each type
    D_cm_kj = diag * gR * (m_k * m_Sigma - m_k_prod.sum(1) * m_Sigma**2 / N + 2. * m_k**2 * m_Sigma**2 * (1. - m_k.sum() * m_Sigma / N))
    D_cm_kj += diag * mR * (1. - m_k.sum() * m_Sigma / N + 2. * (m_k * m_Sigma - m_k_prod.sum(1) * m_Sigma**2 / N))
    D_cm_kj += diag * dR * (m_k * m_Sigma - 2. * m_k**2 * m_Sigma**2) / N
    
    # Co-moment for types k and l
    for k in range(n_types):
        for l in range(n_types):
            if k != l:
            
                D_cm_kj[k,l] += (gR[k] + gR[l]) * m_k_prod[k,l] * m_Sigma**2 * (1. - m_k.sum() * m_Sigma / N)
                D_cm_kj[k,l] += mR[k] * (m_k[l] * m_Sigma - m_k_prod[l,:].sum() * m_Sigma**2 / N) + mR[l] * (m_k[k] * m_Sigma - m_k_prod[k,:].sum() * m_Sigma**2 / N)
                D_cm_kj[k,l] -= (dR[k] + dR[l]) * m_k_prod[k,l] * m_Sigma**2 / N
    
    # Second moment for each type (continuation)
    D_cm_kj -= diag * 2. * m_k**2 * m_Sigma * D_m_Sigma
    
    # Co-moment for types k and j (continuation)
    D_cm_kj -= off_diag * 2. * m_k_prod * m_Sigma * D_m_Sigma
    
    D_cm_kj /= m_Sigma**2

    return np.hstack((D_m_k, D_cm_kj.reshape(1, n_types**2)[0], D_m_Sigma))


# Stopping event for diverging numerical solutions
def stopping_event_abs_abund(t, y, gR, dR, mR, N):
    return all(y[:n_types] - upper_threshold < 0)
stopping_event_abs_abund.terminal = True

# Numerical solver for a set of parameters (absolute abundance)
def model_abs_abund_closure_from_3rd(parameters):
    
    # Growth rates
    gR = np.array([parameters['gR_%i'%i] for i in range(n_types)])
    # Death rates
    dR = np.array([parameters['dR_%i'%i] for i in range(n_types)])
    # Immigration rates
    mR = np.array([parameters['mR_%i'%i] for i in range(n_types)])
    # Carrying capacity
    N = parameters['N']

    # Check that all parameters values are meaningful
    if sum((gR<0) + (dR<0) + (mR<0)) + (N<0) == 0:
        
        # Solve the model numerically
        predicted_moments = integrate.solve_ivp(eq_moments_abs_abund_closure_from_3rd, [0, t_simulated], init_abs_abund, args = (gR, dR, mR, N), t_eval = sampling_times, method = 'LSODA', events=stopping_event_abs_abund)

        # Check if a solution was found
        if predicted_moments.success and predicted_moments.y.shape == (n_types**2 + n_types, len(sampling_times)):

            # Solution
            predicted_moments = predicted_moments.y.T

            return {"moments": predicted_moments}
        
        # Penalize the parameter set if a solution was not found
        else:

            return {"moments": np.inf * np.ones((len(sampling_times), n_types**2 + n_types))}
    
    # Penalize non-meaninful parameter sets
    else:

        return {"moments": np.inf * np.ones((len(sampling_times), n_types**2 + n_types))}
    

# Numerical solver for a set of parameters (absolute abundance)
def model_abs_abund_closure_from_2nd(parameters):
    
    # Growth rates
    gR = np.array([parameters['gR_%i'%i] for i in range(n_types)])
    # Death rates
    dR = np.array([parameters['dR_%i'%i] for i in range(n_types)])
    # Immigration rates
    mR = np.array([parameters['mR_%i'%i] for i in range(n_types)])
    # Carrying capacity
    N = parameters['N']

    # Check that all parameters values are meaningful
    if sum((gR<0) + (dR<0) + (mR<0)) + (N<0) == 0:
        
        # Solve the model numerically
        predicted_moments = integrate.solve_ivp(eq_moments_abs_abund_closure_from_2nd, [0, t_simulated], init_abs_abund, args = (gR, dR, mR, N), t_eval = sampling_times, method = 'LSODA', events=stopping_event_abs_abund)

        # Check if a solution was found
        if predicted_moments.success and predicted_moments.y.shape == (n_types**2 + n_types, len(sampling_times)):

            # Solution
            predicted_moments = predicted_moments.y.T

            return {"moments": predicted_moments}
        
        # Penalize the parameter set if a solution was not found
        else:

            return {"moments": np.inf * np.ones((len(sampling_times), n_types**2 + n_types))}
    
    # Penalize non-meaninful parameter sets
    else:

        return {"moments": np.inf * np.ones((len(sampling_times), n_types**2 + n_types))}

    
# Stopping event for diverging numerical solutions
def stopping_event_rel_abund(t, y, gR, dR, mR, N):
    return y[-1] - 1E6 < 0#True
stopping_event_rel_abund.terminal = True

# Numerical solver for a set of parameters (relative abundance)
def model_rel_abund_closure_from_3rd(parameters):

    # Growth rates
    gR = np.array([parameters['gR_%i'%i] for i in range(n_types)])
    # Death rates
    dR = np.array([parameters['dR_%i'%i] for i in range(n_types)])
    # Immigration rates
    mR = np.array([parameters['mR_%i'%i] for i in range(n_types)])
    # Carrying capacity
    N = parameters['N']
    # First moment of initial abundance
    m_Sigma = parameters['mSigma[0]']

    # Check that all parameters values are meaningful
    if sum((gR<0) + (dR<0) + (mR<0)) + (N<0) + (m_Sigma<0) == 0:

        # Solve the model numerically
        predicted_moments = integrate.solve_ivp(eq_moments_rel_abund_closure_from_3rd, [0, t_simulated], np.hstack((init_rel_abund, m_Sigma)), args = (gR, dR, mR, N), t_eval = sampling_times, method = 'LSODA', events=stopping_event_rel_abund)
        
        # Check if a solution was found
        if predicted_moments.success and predicted_moments.y.shape == (n_types**2 + n_types + 1, len(sampling_times)):
            
            # Solution
            predicted_moments = predicted_moments.y.T

            return {"moments": predicted_moments}

        # Penalize the parameter set if a solution was not found
        else:

            return {"moments": np.inf * np.ones((len(sampling_times), n_types**2 + n_types + 1))}
    
    # Penalize non-meaninful parameter sets
    else:

        return {"moments": np.inf * np.ones((len(sampling_times), n_types**2 + n_types + 1))}
    

# Numerical solver for a set of parameters (relative abundance)
def model_rel_abund_closure_from_2nd(parameters):

    # Growth rates
    gR = np.array([parameters['gR_%i'%i] for i in range(n_types)])
    # Death rates
    dR = np.array([parameters['dR_%i'%i] for i in range(n_types)])
    # Immigration rates
    mR = np.array([parameters['mR_%i'%i] for i in range(n_types)])
    # Carrying capacity
    N = parameters['N']
    # First moment of initial abundance
    m_Sigma = parameters['mSigma[0]']

    # Check that all parameters values are meaningful
    if sum((gR<0) + (dR<0) + (mR<0)) + (N<0) + (m_Sigma<0) == 0:

        # Solve the model numerically
        predicted_moments = integrate.solve_ivp(eq_moments_rel_abund_closure_from_2nd, [0, t_simulated], np.hstack((init_rel_abund, m_Sigma)), args = (gR, dR, mR, N), t_eval = sampling_times, method = 'LSODA', events=stopping_event_rel_abund)
        
        # Check if a solution was found
        if predicted_moments.success and predicted_moments.y.shape == (n_types**2 + n_types + 1, len(sampling_times)):
            
            # Solution
            predicted_moments = predicted_moments.y.T

            return {"moments": predicted_moments}

        # Penalize the parameter set if a solution was not found
        else:

            return {"moments": np.inf * np.ones((len(sampling_times), n_types**2 + n_types + 1))}
    
    # Penalize non-meaninful parameter sets
    else:

        return {"moments": np.inf * np.ones((len(sampling_times), n_types**2 + n_types + 1))}