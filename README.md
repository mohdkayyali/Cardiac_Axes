# Cardiac_Axes

This repository contains Python code for analysing cardiac anatomical and electrical axes. It focuses on extracting axes based on multiple definitions and converting vectors into polar coordinates. The repository is organised into distinct modules, each serving a specific purpose.

## Directory Structure

- **Anatomical_axis/**  
  Contains code to extract anatomical axes based on five different definitions. These calculations are performed using surface meshes with labelled cardiac structures.

- **Electrical_axis/**  
  Includes code to process read ECG signals from xml files, convert them into vectorcardiograms (VCG), and extract various ECG parameters from the xml file. The file also provides functions to compute electrical axes using different definitions.

- **Processing_axes**  
  Scripts to extract anato-electrical angular separation in 3D as well as anatomical and electrical axis orientation metrics in 2D and in the spherical coordinate system. 

- **axis_separation_metrics**  
  MATLAB scripts for computing the error metrics used to identify the optimal cardiac axes definitions. 

## Features

1. **Anatomical Axis Extraction**  
   - Processes 3D surface meshes to determine anatomical axes -- from vtk files.  
   - Supports five predefined methods for axis determination.
   - Computes spherecity index of the left ventricle. 

2. **Electrical Axis Analysis**
   - Extracts ECG signals and parameters from the xml files.  
   - Converts standard ECG signals to VCG.   
   - Computes electrical axes using  distinct definitions.  

3. **Processing_axes**  
   - Transforms 3D vector data into spherical coordinates for further analysis.
   - Computes 3D angular sepration between anatomical & electrical axes
   - Computes 2D, planar angles for individual axes.
