#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 15:15:29 2024

@author: mka23
"""

#######################################################################################
##This file extracts anatomical axes using 5 methods and saves them into a .csv file###
#######################################################################################
# Extracting ES & ED 


import os
import pandas as pd
from anatomical_axis_computation_functions import axis_by_PCA,axis_by_labels, spherecity_index
import pyvista as pv
import concurrent.futures
import re

##Identify files with meshes
source_dir = '/netapp/cme_digital_twins/UKBB_88878/Data/'  
files = os.listdir(source_dir)
Mesh_files = [file for file in files if file.startswith('UKBB')]

def import_mesh_main(datadir, timeframe):
    mesh_TV = pv.read(datadir + f'/Mesh_tricuspid_valve_timeframe{timeframe}.vtk')
    mesh_PV = pv.read(datadir + f'/Mesh_pulmonary_valve_timeframe{timeframe}.vtk')
    mesh_AV = pv.read(datadir + f'/Mesh_aorta_valve_timeframe{timeframe}.vtk')
    mesh_MV = pv.read(datadir + f'/Mesh_mitral_valve_timeframe{timeframe}.vtk')
    mesh_RV_FW = pv.read(datadir + f'/Mesh_RV_FW_timeframe{timeframe}.vtk')
    mesh_RV_septum = pv.read(datadir + f'/Mesh_RV_septum_timeframe{timeframe}.vtk')
    mesh_LV_endo = pv.read(datadir + f'/Mesh_LV_endo_timeframe{timeframe}.vtk')
    mesh_epi = pv.read(datadir + f'/Mesh_epi_timeframe{timeframe}.vtk')
    
    return {
        'TV': mesh_TV,
        'PV': mesh_PV,
        'AV': mesh_AV,
        'MV': mesh_MV,
        'RV_FW': mesh_RV_FW,
        'RV_septum': mesh_RV_septum,
        'LV_endo': mesh_LV_endo,
        'epi': mesh_epi
    }

def find_timeframes(directory):
    filenames = os.listdir(directory)
    timeframes = set()
    pattern = re.compile(r'timeframe(\d{3})\.vtk')

    for filename in filenames:
        match = pattern.search(filename)
        if match:
            timeframes.add(match.group(1))

    sorted_timeframes = sorted(timeframes)
    return sorted_timeframes


# Function to process a single patient file
def process_patient_file(patient_file):
    patient_id = patient_file.split('_')[2]
    source_file_path = os.path.join(source_dir, patient_file)
    path = os.path.join(source_file_path, 'Instance_2')
    
    try:
        patient_folder = os.path.join(path, 'Mesh_Outputs')
        timeframes = find_timeframes(patient_folder)

        if len(timeframes) == 0:
            return None, None  # Skip if no timeframes found

        mesh_dict_ED = import_mesh_main(patient_folder, timeframes[0])
        mesh_dict_ES = import_mesh_main(patient_folder, timeframes[-1]) if len(timeframes) > 1 else None

        def calculate_metrics(mesh_dict, timeframe):
            # Calculate axis by PCA method for each mesh part
            axis_PCA_LV = axis_by_PCA(mesh_dict, 'LV_endo')
            axis_PCA_All = axis_by_PCA(mesh_dict, 'Whole')

            # Calculate axis by labels for specified mesh parts
            base_key_1 = ['MV']
            base_key_2 = ['MV', 'AV']
            base_key_3 = ['MV', 'AV', 'TV', 'PV']
            apex_type = 'endo'
            axis_labels_MV, apex_MV, base_MV = axis_by_labels(mesh_dict, base_key_1, apex_type)
            axis_labels_MV_AV, apex, base_MV_AV = axis_by_labels(mesh_dict, base_key_2, apex_type)
            axis_labels_AllV, apex_AllV, base_AllV = axis_by_labels(mesh_dict, base_key_3, apex_type)

            Spherecity = spherecity_index(mesh_dict, apex_MV, base_MV)

            LV_Volume = mesh_dict['LV_endo'].volume
            Epi_Volume = mesh_dict['epi'].volume

            # Return the results as a dictionary
            return {'eid': patient_id,
                    'Timeframe': timeframe,
                    'Axis_PCA_LV': axis_PCA_LV, 'Axis_PCA_All': axis_PCA_All,
                    'Axis_Labels_MV': axis_labels_MV, 'Axis_Labels_MV_AV': axis_labels_MV_AV,
                    'Axis_Labels_AllV': axis_labels_AllV,
                    'Spherecity': Spherecity,
                    'Volume_LV': LV_Volume,
                    'Volume_Epi': Epi_Volume}

        result_ED = calculate_metrics(mesh_dict_ED, timeframes[0])
        result_ES = calculate_metrics(mesh_dict_ES, timeframes[-1]) if mesh_dict_ES else None

        return result_ED, result_ES

    except FileNotFoundError:
        print(f"Mesh parts not found for patient {patient_id}")
    except Exception as e:
        print(f"An unexpected error occurred for patient {patient_id}: {e}")
    return None, None



#%% Main processing loop with parallel execution
results_ED = []
results_ES = []

with concurrent.futures.ProcessPoolExecutor() as executor:
    futures = [executor.submit(process_patient_file, patient_file) for patient_file in Mesh_files]
    for future in concurrent.futures.as_completed(futures):
        result_ED, result_ES = future.result()
        if result_ED:
            results_ED.append(result_ED)
        if result_ES:
            results_ES.append(result_ES)

# Convert the lists of results to DataFrames
result_df_ED = pd.DataFrame(results_ED)
result_df_ES = pd.DataFrame(results_ES)

result_df_ED.to_csv('anatomical_axis_results_full_ED.csv', index=False)  # Save without row indices
result_df_ES.to_csv('anatomical_axis_results_full_ES.csv', index=False)  # Save without row indices

