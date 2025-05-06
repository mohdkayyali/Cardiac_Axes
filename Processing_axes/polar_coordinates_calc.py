#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri July  19 11:19:55 2024

@author: mka23
"""

import math
import numpy as np
import pandas as pd
from vector_meta_data import vector_data_all, vector_data_healthy
from concurrent.futures import ProcessPoolExecutor

def cartesian_to_polar(vector_3D):
    x, z, y = vector_3D  # z and y are switched here to allow frontal plane (ZX) angle == theta
    r = math.sqrt(x**2 + y**2 + z**2)
    theta = math.atan2(y, x)
    phi = math.acos(z / r)
    return r, math.degrees(theta), math.degrees(phi)

def process_chunk(chunk, inverted_anat=True):
    polar_coordinates_chunk = pd.DataFrame(columns=['eid', 'theta_A', 'phi_A', 'theta_E', 'phi_E'])
    
    for index, row in chunk.iterrows():
        eid = str(row['eid'])
        vector_A = row['AM5']
        vector_E = row['EM1']
        if inverted_anat:
            vector_A = -vector_A
        vector_A /= np.linalg.norm(vector_A)
        vector_E /= np.linalg.norm(vector_E)
        
        _, theta_A, phi_A = cartesian_to_polar(vector_A)
        _, theta_E, phi_E = cartesian_to_polar(vector_E)
        
        polar_coordinates_chunk.loc[index] = [eid, theta_A, phi_A, theta_E, phi_E]
        
    return polar_coordinates_chunk

def polar_coords_all_vectors_parallel(vector_data, inverted_anat=True, num_workers=4):
    chunk_size = int(math.ceil(len(vector_data) / num_workers))
    chunks = [vector_data.iloc[i:i + chunk_size] for i in range(0, len(vector_data), chunk_size)]
    
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        results = executor.map(process_chunk, chunks, [inverted_anat] * num_workers)
    
    polar_coordinates = pd.concat(results).reset_index(drop=True)
    
    return polar_coordinates

# Run parallel processing for vector_data_healthy
polar_coordinates_healthy_inv = polar_coords_all_vectors_parallel(vector_data_healthy, inverted_anat=True)
polar_coordinates_healthy = polar_coords_all_vectors_parallel(vector_data_healthy, inverted_anat=False)

polar_coordinates_healthy['dtheta'] = polar_coordinates_healthy['theta_A'] - polar_coordinates_healthy['theta_E']
polar_coordinates_healthy['dphi'] = polar_coordinates_healthy['phi_A'] - polar_coordinates_healthy['phi_E']
polar_coordinates_healthy_inv['dtheta'] = polar_coordinates_healthy_inv['theta_A'] - polar_coordinates_healthy_inv['theta_E']
polar_coordinates_healthy_inv['dphi'] = polar_coordinates_healthy_inv['phi_A'] - polar_coordinates_healthy_inv['phi_E']

polar_coordinates_healthy.loc[polar_coordinates_healthy['dtheta'] < -180, 'dtheta'] += 360
polar_coordinates_healthy.loc[polar_coordinates_healthy['dphi'] < -180, 'dtheta'] += 360
polar_coordinates_healthy_inv.loc[polar_coordinates_healthy_inv['dtheta'] < -180, 'dtheta'] += 360
polar_coordinates_healthy_inv.loc[polar_coordinates_healthy_inv['dphi'] < -180, 'dtheta'] += 360

# Run parallel processing for vector_data_all
polar_coordinates_all_inv = polar_coords_all_vectors_parallel(vector_data_all, inverted_anat=True)
polar_coordinates_all = polar_coords_all_vectors_parallel(vector_data_all, inverted_anat=False)

polar_coordinates_all['dtheta'] = polar_coordinates_all['theta_A'] - polar_coordinates_all['theta_E']
polar_coordinates_all['dphi'] = polar_coordinates_all['phi_A'] - polar_coordinates_all['phi_E']
polar_coordinates_all_inv['dtheta'] = polar_coordinates_all_inv['theta_A'] - polar_coordinates_all_inv['theta_E']
polar_coordinates_all_inv['dphi'] = polar_coordinates_all_inv['phi_A'] - polar_coordinates_all_inv['phi_E']

polar_coordinates_all.loc[polar_coordinates_all['dtheta'] < -180, 'dtheta'] += 360
polar_coordinates_all.loc[polar_coordinates_all['dphi'] < -180, 'dtheta'] += 360
polar_coordinates_all_inv.loc[polar_coordinates_all_inv['dtheta'] < -180, 'dtheta'] += 360
polar_coordinates_all_inv.loc[polar_coordinates_all_inv['dphi'] < -180, 'dtheta'] += 360

polar_coordinates_healthy.to_csv('polar_coordinates_results/polar_coordinates_healthy.csv', index=False)
polar_coordinates_healthy_inv.to_csv('polar_coordinates_results/polar_coordinates_healthy_inv.csv', index=False)
polar_coordinates_all.to_csv('polar_coordinates_results/polar_coordinates_all.csv', index=False)
polar_coordinates_all_inv.to_csv('polar_coordinates_results/polar_coordinates_all_inv.csv', index=False)
