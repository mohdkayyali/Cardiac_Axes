#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 22:19:54 2024
Modified on May 1 2025 to include both 3D and full 5-method planar outputs

@author: mka23 (modified)
"""

import math
import numpy as np
import pandas as pd
from vector_meta_data import vector_data_all, vector_data_healthy
from concurrent.futures import ProcessPoolExecutor

# ─── Helpers ──────────────────────────────────────────────────────────────────

def calculate_angle_and_cosine(vector_A, vector_E):
    vector_A = vector_A / np.linalg.norm(vector_A)
    vector_E = vector_E / np.linalg.norm(vector_E)
    dot = np.clip(np.dot(vector_A, vector_E), -1.0, 1.0)
    angle = math.degrees(math.acos(dot))
    return angle, dot

def planar_angles(vec):
    x, y, z = vec
    return (
        math.degrees(math.atan2(z, x)),  # α (X–Z)
        math.degrees(math.atan2(z, y)),  # β (Y–Z)
        math.degrees(math.atan2(y, x)),  # γ (X–Y)
    )

def wrap_diff(arr):
    arr = np.where(arr < -180, arr + 360, arr)
    arr = np.where(arr >  180, arr - 360, arr)
    return arr

# ─── 3D chunk processor ─────────────────────────────────────────────────────

def process_chunk_3d(args):
    chunk, inverted_anat = args
    rows = []
    for _, row in chunk.iterrows():
        eid = str(row['eid'])
        A = row['AM5'].astype(float)
        E = row['EM1'].astype(float)
        if inverted_anat:
            A = -A
        angle3d, cos3d = calculate_angle_and_cosine(A, E)
        rows.append({
            'eid':          eid,
            'AE_3D_angle':  angle3d,
            'AE_3D_cosine': cos3d
        })
    return pd.DataFrame(rows)

# ─── Planar chunk processor ─────────────────────────────────────────────────

def process_chunk_planar(args):
    chunk, inverted_anat = args
    rows = []
    for _, row in chunk.iterrows():
        eid = str(row['eid'])

        # grab and invert if needed
        E_vecs = [row[f'EM{i}'].astype(float) for i in range(1,6)]
        A_vecs = [(-row[f'AM{i}'].astype(float) if inverted_anat 
                   else row[f'AM{i}'].astype(float)) for i in range(1,6)]

        # per-method α/β/γ
        alpha_E = []; beta_E = []; gamma_E = []
        alpha_A = []; beta_A = []; gamma_A = []
        for E, A in zip(E_vecs, A_vecs):
            aE, bE, gE = planar_angles(E)
            aA, bA, gA = planar_angles(A)
            alpha_E.append(aE); beta_E.append(bE); gamma_E.append(gE)
            alpha_A.append(aA); beta_A.append(bA); gamma_A.append(gA)

        # 25 deltas for each angle type (A minus E), wrapped
        flat_da = []; flat_db = []; flat_dg = []
        for i in range(5):
            for j in range(5):
                flat_da.append(wrap_diff(np.array([alpha_A[j] - alpha_E[i]]))[0])
                flat_db.append(wrap_diff(np.array([beta_A[j]  - beta_E[i]]))[0])
                flat_dg.append(wrap_diff(np.array([gamma_A[j] - gamma_E[i]]))[0])

        # assemble row
        row_dict = {'eid': eid}
        # alpha_A_1…5, beta_A_1…5, gamma_A_1…5
        for i in range(5):
            row_dict[f'alpha_A_{i+1}'] = alpha_A[i]
            row_dict[f'beta_A_{i+1}']  = beta_A[i]
            row_dict[f'gamma_A_{i+1}'] = gamma_A[i]
        # alpha_E_1…5, beta_E_1…5, gamma_E_1…5
        for i in range(5):
            row_dict[f'alpha_E_{i+1}'] = alpha_E[i]
            row_dict[f'beta_E_{i+1}']  = beta_E[i]
            row_dict[f'gamma_E_{i+1}'] = gamma_E[i]
        # delta_alpha_1…25
        for k, v in enumerate(flat_da, start=1):
            row_dict[f'delta_alpha_{k}'] = v
        # delta_beta_1…25
        for k, v in enumerate(flat_db, start=1):
            row_dict[f'delta_beta_{k}'] = v
        # delta_gamma_1…25
        for k, v in enumerate(flat_dg, start=1):
            row_dict[f'delta_gamma_{k}'] = v

        rows.append(row_dict)

    return pd.DataFrame(rows)

# ─── Parallel wrappers ────────────────────────────────────────────────────────

def ae_3d_all_vectors_parallel(vector_data, inverted_anat=True, num_workers=4):
    chunk_size = int(math.ceil(len(vector_data) / num_workers))
    tasks = [
        (vector_data.iloc[i:i+chunk_size], inverted_anat)
        for i in range(0, len(vector_data), chunk_size)
    ]
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        dfs = executor.map(process_chunk_3d, tasks)
    return pd.concat(dfs).reset_index(drop=True)

def ae_planar_all_vectors_parallel(vector_data, inverted_anat=True, num_workers=4):
    chunk_size = int(math.ceil(len(vector_data) / num_workers))
    tasks = [
        (vector_data.iloc[i:i+chunk_size], inverted_anat)
        for i in range(0, len(vector_data), chunk_size)
    ]
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        dfs = executor.map(process_chunk_planar, tasks)
    return pd.concat(dfs).reset_index(drop=True)

# ─── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    # Healthy
    AE_3D_healthy     = ae_3d_all_vectors_parallel(vector_data_healthy, inverted_anat=False, num_workers=4)
    AE_planar_healthy = ae_planar_all_vectors_parallel(vector_data_healthy, inverted_anat=False, num_workers=4)
    AE_3D_healthy.to_csv('AE_3D_healthy.csv',     index=False)
    AE_planar_healthy.to_csv('AE_2D_healthy.csv', index=False)

    # All
    AE_3D_all     = ae_3d_all_vectors_parallel(vector_data_all, inverted_anat=False, num_workers=4)
    AE_planar_all = ae_planar_all_vectors_parallel(vector_data_all, inverted_anat=False, num_workers=4)
    AE_3D_all.to_csv('AE_3D_all.csv',     index=False)
    AE_planar_all.to_csv('AE_2D_all.csv', index=False)

    print("Done: 3D (angle+cosine) and planar (5 methods + 25 deltas) CSVs written.")
