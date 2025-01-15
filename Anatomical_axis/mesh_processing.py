# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 12:48:41 2023

@author: mkayy
"""

from sklearn.preprocessing import StandardScaler
from functions import points_from_faces
import pyvista as pv
import numpy as np
def import_mesh_2(datadir):
    mesh_TV = pv.read(datadir + '/Mesh_tricuspid_valve_timeframe001.vtk')
    mesh_PV = pv.read(datadir + '/Mesh_pulmonary_valve_timeframe001.vtk')
    mesh_AV = pv.read(datadir + '/Mesh_aorta_valve_timeframe001.vtk')
    mesh_MV = pv.read(datadir + '/Mesh_mitral_valve_timeframe001.vtk')
    mesh_RV_FW = pv.read(datadir + '/Mesh_RV_FW_timeframe001.vtk')
    mesh_RV_septum = pv.read(datadir + '/Mesh_RV_septum_timeframe001.vtk')
    mesh_LV_endo = pv.read(datadir + '/Mesh_LV_endo_timeframe001.vtk')
    mesh_epi = pv.read(datadir + '/Mesh_epi_timeframe001.vtk')
    
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

def import_mesh(patient_2964_ID, datadir):
    datadir = datadir + f'2964_{patient_2964_ID}/Mesh_Outputs/'
    mesh_TV = pv.read(datadir + 'Mesh_tricuspid_valve_timeframe001.vtk')
    mesh_PV = pv.read(datadir + 'Mesh_pulmonary_valve_timeframe001.vtk')
    mesh_AV = pv.read(datadir + 'Mesh_aorta_valve_timeframe001.vtk')
    mesh_MV = pv.read(datadir + 'Mesh_mitral_valve_timeframe001.vtk')
    mesh_RV_FW = pv.read(datadir + 'Mesh_RV_FW_timeframe001.vtk')
    mesh_RV_septum = pv.read(datadir + 'Mesh_RV_septum_timeframe001.vtk')
    mesh_LV_endo = pv.read(datadir + 'Mesh_LV_endo_timeframe001.vtk')
    mesh_epi = pv.read(datadir + 'Mesh_epi_timeframe001.vtk')
    
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

def axis_by_PCA(mesh_dict, mesh_key):
    if mesh_key == 'Whole':
        points_mesh = mesh_dict['epi'].points #doesn't matter which part, all have same points (diff faces)
    else:
        # Extract points from faces
        specific_points = points_from_faces(mesh_dict[mesh_key].faces)
        points_mesh = np.unique(np.array([mesh_dict[mesh_key].points[i] for i in specific_points]), axis=0)
    # Center the data
    scaler = StandardScaler()
    standardized_data = scaler.fit_transform(points_mesh)

    # Perform PCA
    covariance_matrix = np.cov(standardized_data, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance_matrix)

    # Identify the long axis
    axis = eigenvectors[:, -1]

    return axis

def axis_by_labels(mesh_dict, base_mesh_key, apex_type):
    
    centroids_list = []
    for i in range(len(base_mesh_key)):   
        cell_centers = mesh_dict[base_mesh_key[i]].cell_centers()
        centroids_list.append(np.mean(cell_centers.points, axis=0))

    centroid = np.array(centroids_list)       
    base = np.mean(centroid, axis=0)

    if apex_type == 'endo':
        specific_points = points_from_faces(mesh_dict['LV_endo'].faces)
        points_specific_mesh = np.unique(np.array([mesh_dict['LV_endo'].points[i] for i in specific_points]), axis=0)

        # Calculate the distances from the base to all points in the mesh
        distances_to_apex = np.linalg.norm(points_specific_mesh - base, axis=1)

        # Get the coordinates of the point with the maximum distance
        max_distance_idx = np.argmax(distances_to_apex)
        apex = points_specific_mesh[max_distance_idx]
        axis = base-apex
        return axis, apex, base
    
        
    elif apex_type == 'epi':
        specific_points = points_from_faces(mesh_dict['epi'].faces)
        points_specific_mesh = np.unique(np.array([mesh_dict['epi'].points[i] for i in specific_points]), axis=0)

        # Calculate the distances from the base to all points in the mesh
        distances_to_apex = np.linalg.norm(points_specific_mesh - base, axis=1)

        # Get the coordinates of the point with the maximum distance
        max_distance_idx = np.argmax(distances_to_apex)
        apex = points_specific_mesh[max_distance_idx]
        axis = base-apex
        return axis, apex, base
    
    else:
        raise ValueError("Invalid 'apex_type'. Use 'endo' or 'epi'.")


def spherecity_index(mesh_dict, apex, base):
    """

    Parameters
    ----------
    mesh_dict : Dictionary
        the imported mesh from import_mesh()
    apex : Array
        apex as computed from axis_by_labels function
    base : Array
        base as computed from axis_by_labels function

    Returns
    -------
    spherecity : Float
        Value of spherecity index = (short axis length) / (long axis length)

    """
    LA = np.linalg.norm(base - apex)
    midpoint = (apex + base) / 2

    arrow_direction = (base - apex) / np.linalg.norm(base - apex)

    perpendicular_vector = np.cross(arrow_direction, [1, 0, 0])  # Assumes arrow not aligned with x-axis
    # st = midpoint - 50 * perpendicular_vector
    # end = midpoint + 50 * perpendicular_vector
    # perpendicular_line = pv.Line(midpoint - 50 * perpendicular_vector, midpoint + 50 * perpendicular_vector)

    # Compute multiple perp vectors
    for angle_degrees in range(0, 180, 5):
        angle_radians = np.radians(angle_degrees)

        # Rotate the perpendicular vector using a rotation matrix (Rodrigues' Rotation Formula)
        rotation_matrix = np.array([
            [np.cos(angle_radians) + arrow_direction[0]**2 * (1 - np.cos(angle_radians)),
             arrow_direction[0] * arrow_direction[1] * (1 - np.cos(angle_radians)) - arrow_direction[2] * np.sin(angle_radians),
             arrow_direction[0] * arrow_direction[2] * (1 - np.cos(angle_radians)) + arrow_direction[1] * np.sin(angle_radians)],
            [arrow_direction[1] * arrow_direction[0] * (1 - np.cos(angle_radians)) + arrow_direction[2] * np.sin(angle_radians),
             np.cos(angle_radians) + arrow_direction[1]**2 * (1 - np.cos(angle_radians)),
             arrow_direction[1] * arrow_direction[2] * (1 - np.cos(angle_radians)) - arrow_direction[0] * np.sin(angle_radians)],
            [arrow_direction[2] * arrow_direction[0] * (1 - np.cos(angle_radians)) - arrow_direction[1] * np.sin(angle_radians),
             arrow_direction[2] * arrow_direction[1] * (1 - np.cos(angle_radians)) + arrow_direction[0] * np.sin(angle_radians),
             np.cos(angle_radians) + arrow_direction[2]**2 * (1 - np.cos(angle_radians))]
        ])

        perpendicular_vector_rotated = np.dot(rotation_matrix, perpendicular_vector)

        intersecting_cells = mesh_dict['LV_endo'].find_cells_intersecting_line(midpoint - 100 * perpendicular_vector_rotated, midpoint + 100 * perpendicular_vector_rotated)
        #int_cells = mesh_dict['LV_endo'].extract_cells(intersecting_cells)
        #corresponding_points = int_cells.points
        n_cells = len(intersecting_cells)
        centroids = np.zeros((n_cells,3))
        for i in range(0,n_cells):
            centroids[i] = np.mean(mesh_dict['LV_endo'].extract_cells(intersecting_cells[i]).points, axis=0)
        if n_cells >2:
            print('There is more than 1 cells intersecting with the line at each side for:' + f'{angle_degrees}')
                   
        current_distance = sum(np.linalg.norm(centroids - midpoint, axis=1))
        if angle_degrees == 0:
            min_distance = current_distance
        if current_distance < min_distance:
            min_distance = current_distance
            # angle_at_min_distance = angle_degrees
            # SA_centroids = centroids

    SA = min_distance
    
    spherecity = SA/LA
    
    return spherecity