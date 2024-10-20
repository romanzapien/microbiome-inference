#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@project: microbiome inference (lotka-volterra model - source code)
@author: Roman Zapien-Campos - 2023
"""

# Import packages
import numpy as np
from scipy import integrate
from par_LV_numerics import *


# Indexes of elements in a matrix of size n_types
diag = np.eye(n_types)
off_diag = 1 - diag


## Moment Closure - 3rd moments & one 2nd moment (WORKING - AFTER TRYING MANY SOLUTIONS)
def eq_moments_abs_abund_closure_from_3rd(t, y, gR, I_p, I_n): 
    
    # Extract different moments
    m_k = y[:n_types]
    cm_kj = y[n_types:].reshape(n_types,n_types)
    m_kk = np.diag(cm_kj)
    m_k_prod = np.outer(m_k, m_k)
    
    # First moment for each type
    D_m_k = gR * m_k
    D_m_k += ((I_p - I_n) * m_k_prod).sum(1)
    
    # Second moment for each type
    D_cm_kj = diag * gR * (m_k + 2. * m_kk)
    D_cm_kj += diag * ((I_p + I_n) * cm_kj).sum(1)
    D_cm_kj += diag * 2. * ((I_p - I_n) * m_k).sum(1) * m_kk
    
    # Co-moment for types k and l
    for k in range(n_types):
        for l in range(n_types):
            if k != l:
            
                D_cm_kj[k,l] += (gR[k] + gR[l]) * cm_kj[k,l]
                D_cm_kj[k,l] += ((I_p[k,:] - I_n[k,:] + I_p[l,:] - I_n[l,:]) * cm_kj[k,l] * m_k).sum()

    return np.hstack((D_m_k, D_cm_kj.reshape(1, n_types**2)[0]))


## Moment Closure - 3rd moments & one 2nd moment (WORKING - AFTER TRYING MANY SOLUTIONS)
def eq_moments_rel_abund_closure_from_3rd(t, y, gR, I_p, I_n):

    # Unpack different moments
    m_k = y[:n_types]
    cm_kj = y[n_types:-1].reshape(n_types,n_types)
    m_kk = np.diag(cm_kj)
    m_Sigma = y[-1]
    m_k_prod = np.outer(m_k, m_k)

    # First moment for each type
    D_m_k = gR * m_k * m_Sigma
    D_m_k += ((I_p - I_n) * m_k_prod * m_Sigma**2).sum(1)

    # First moment of scaling factor
    D_m_Sigma = sum(D_m_k)

    # First moment for each type (continuation)
    D_m_k -= m_k * D_m_Sigma
    
    D_m_k /= m_Sigma

    # Second moment for each type
    D_cm_kj = diag * gR * (m_k * m_Sigma + 2. * m_kk * m_Sigma**2)
    D_cm_kj += diag * ((I_p + I_n) * cm_kj * m_Sigma**2).sum(1)
    D_cm_kj += diag * 2. * ((I_p - I_n) * m_k * m_Sigma).sum(1) * m_kk * m_Sigma**2

    # Co-moment for types k and j
    for k in range(n_types):
        for l in range(n_types):
            if k != l:

                D_cm_kj[k,l] += (gR[k] + gR[l]) * cm_kj[k,l] * m_Sigma**2
                D_cm_kj[k,l] += ((I_p[k,:] - I_n[k,:] + I_p[l,:] - I_n[l,:]) * cm_kj[k,l] * m_k * m_Sigma**3).sum()

    # Second moment for each type (continuation)
    D_cm_kj -= diag * 2. * m_kk * m_Sigma * D_m_Sigma

    # Co-moment for types k and j (continuation)
    D_cm_kj -= off_diag * 2. * cm_kj * m_Sigma * D_m_Sigma
    
    D_cm_kj /= m_Sigma**2

    return np.hstack((D_m_k, D_cm_kj.reshape(1, n_types**2)[0], D_m_Sigma))


## Moment Closure - 2nd & 3rd moments (WORKING)
def eq_moments_abs_abund_closure_from_2nd(t, y, gR, I_p, I_n): 

   # Unpack different moments
    m_k = y[:n_types]
    m_k_prod = np.outer(m_k, m_k)
    
    # First moment for each type
    D_m_k = gR * m_k
    D_m_k += ((I_p - I_n) * m_k_prod).sum(1)
    
    # Second moment for each type
    D_cm_kj = diag * gR * (m_k + 2. * m_k**2)
    D_cm_kj += diag * ((I_p + I_n) * m_k_prod).sum(1)
    D_cm_kj += diag * 2. * ((I_p - I_n) * m_k_prod).sum(1) * m_k
    
    # Co-moment for types k and l
    for k in range(n_types):
        for l in range(n_types):
            if k != l:
            
                D_cm_kj[k,l] += (gR[k] + gR[l]) * m_k_prod[k,l]
                D_cm_kj[k,l] += ((I_p[k,:] - I_n[k,:] + I_p[l,:] - I_n[l,:]) * m_k_prod[k,l] * m_k).sum()

    return np.hstack((D_m_k, D_cm_kj.reshape(1, n_types**2)[0]))


## Moment Closure - 2nd & 3rd moments (WORKING)
def eq_moments_rel_abund_closure_from_2nd(t, y, gR, I_p, I_n):

   # Unpack different moments
    m_k = y[:n_types]
    m_Sigma = y[-1]
    m_k_prod = np.outer(m_k, m_k)
                
    # First moment for each type
    D_m_k = gR * m_k * m_Sigma
    D_m_k += ((I_p - I_n) * m_k_prod * m_Sigma**2).sum(1)
    
    # First moment of scaling factor
    D_m_Sigma = sum(D_m_k)
    
    # First moment for each type (continuation)
    D_m_k -= m_k * D_m_Sigma
    
    D_m_k /= m_Sigma
    
    # Second moment for each type
    D_cm_kj = diag * gR * (m_k * m_Sigma + 2. * m_k**2 * m_Sigma**2)
    D_cm_kj += diag * ((I_p + I_n) * m_k_prod * m_Sigma**2).sum(1)
    D_cm_kj += diag * 2. * ((I_p - I_n) * m_k * m_Sigma).sum(1) * m_k**2 * m_Sigma**2
    
    # Co-moment for types k and j
    for k in range(n_types):
        for l in range(n_types):
            if k != l:

                D_cm_kj[k,l] += (gR[k] + gR[l]) * m_k_prod[k,l] * m_Sigma**2
                D_cm_kj[k,l] += ((I_p[k,:] - I_n[k,:] + I_p[l,:] - I_n[l,:]) * m_k_prod[k,l] * m_k * m_Sigma**3).sum()
    
    # Second moment for each type (continuation)
    D_cm_kj -= diag * 2. * m_k**2 * m_Sigma * D_m_Sigma
    
    # Co-moment for types k and j (continuation)
    D_cm_kj -= off_diag * 2. * m_k_prod * m_Sigma * D_m_Sigma
    
    D_cm_kj /= m_Sigma**2
        
    return np.hstack((D_m_k, D_cm_kj.reshape(1, n_types**2)[0], D_m_Sigma))


# Stopping event for diverging numerical solutions
def stopping_event_abs_abund(t, y, gR, I_p, I_n):
    return all(y[:n_types] - upper_threshold < 0)
stopping_event_abs_abund.terminal = True

# Numerical solver for a set of parameters (absolute abundance)
def model_abs_abund_closure_from_3rd(parameters):

    # Growth rates
    gR = np.array([parameters['gR_%i'%i] for i in range(n_types)])
    # Interactions
    I = np.array([parameters['I_%i_%i'%(i,j)] for i in range(n_types) for j in range(n_types)]).reshape(n_types,n_types)
    # Intra-specific
    I_intra = I[np.eye(n_types,dtype=bool)]
    # Inter-specific (positive)
    I_p = np.zeros((n_types, n_types))
    I_p[I>0] = abs(I[I>0])
    # Inter-specific (negative)
    I_n = np.zeros((n_types, n_types))
    I_n[I<0] = abs(I[I<0])

    # Check that all parameters values are meaningful
    if sum((gR<0) + (I_intra>=0)) == 0:

        # Solve the model numerically
        predicted_moments = integrate.solve_ivp(eq_moments_abs_abund_closure_from_3rd, [0, t_simulated], init_abs_abund, args = (gR, I_p, I_n), t_eval = sampling_times, method = 'LSODA', events=stopping_event_abs_abund)
        
        # Check if a solution was found
        if predicted_moments.success and predicted_moments.y.shape == (n_types**2 + n_types, t_points):
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
def stopping_event_rel_abund(t, y, gR, I_p, I_n):
    return y[-1] - 1E6 < 0#True
stopping_event_rel_abund.terminal = True

# Numerical solver for a set of parameters (relative abundance)
def model_rel_abund_closure_from_3rd(parameters):

    # Growth rates
    gR = np.array([parameters['gR_%i'%i] for i in range(n_types)])
    # Interactions
    I = np.array([parameters['I_%i_%i'%(i,j)] for i in range(n_types) for j in range(n_types)]).reshape(n_types,n_types)
    # Intra-specific
    I_intra = I[np.eye(n_types,dtype=bool)]
    # Inter-specific (positive)
    I_p = np.zeros((n_types, n_types))
    I_p[I>0] = abs(I[I>0])
    # Inter-specific (negative)
    I_n = np.zeros((n_types, n_types))
    I_n[I<0] = abs(I[I<0])
    # First moment of initial abundance
    m_Sigma = parameters['mSigma[0]']

    # Check that all parameters values are meaningful
    if sum((gR<0) + (I_intra>=0)) + (m_Sigma<0) == 0:
        
        # Solve the model numerically
        predicted_moments = integrate.solve_ivp(eq_moments_rel_abund_closure_from_3rd, [0, t_simulated], np.hstack((init_rel_abund, m_Sigma)), args = (gR, I_p, I_n), t_eval = sampling_times, method = 'LSODA', events=stopping_event_rel_abund)

        # Check if a solution was found
        if predicted_moments.success and predicted_moments.y.shape == (n_types**2 + n_types + 1, t_points):

            # Solution
            predicted_moments = predicted_moments.y.T

            return {"moments": predicted_moments}

        # Penalize the parameter set if a solution was not found
        else:
            
            return {"moments": np.inf * np.ones((len(sampling_times), n_types**2 + n_types + 1))}

    # Penalize non-meaninful parameter sets
    else:

        return {"moments": np.inf * np.ones((len(sampling_times),n_types**2 + n_types + 1))}