#!/usr/bin/env python3
"""
Created on Thu May 16 12:50:27 2019

@author: rishabhstein
"""

import numpy as np

def roughness(dims, f):
    '''Create a random mineral field
        dims - tuple with system dimensions
        f - list of mineral porosity and perturbation'''
    
    data  = np.ones(dims)*f[0]
    data += (2*np.random.random(dims)-1)*np.sqrt(3.)*f[1]
    
    return data

def create_porosity_roughness(dims, porosity = 1, noise = 1):
    '''Create a porosity field with dims, dimensions and
    noise'''
    data  = np.ones(dims)*porosity
    data += (2*np.random.random(dims)-1)*noise
    return data