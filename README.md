# Cardiac_Axes

This repository contains Python code for analysing cardiac anatomical and electrical axes. It focuses on extracting axes based on multiple definitions and converting vectors into polar coordinates. The repository is organised into distinct modules, each serving a specific purpose.

## Directory Structure

- **Anatomical_axis/**  
  Contains code to extract anatomical axes based on five different definitions. These calculations are performed using surface meshes of cardiac structures.

- **Electrical_axis/**  
  Includes code to process electrocardiogram (ECG) signals, convert them into vectorcardiograms (VCG), and extract various ECG parameters. The module also provides methods to determine electrical axes using four different definitions.

- **polar_coordinates_calc.py**  
  Implements functionality to calculate polar coordinates from vector data.

- **polar_coordinates_load.py**  
  A helper script for loading and pre-processing vector data before converting to polar coordinates.

## Features

1. **Anatomical Axis Extraction**  
   - Processes 3D surface meshes to determine anatomical axes.  
   - Supports five predefined methods for axis determination.  

2. **Electrical Axis Analysis**  
   - Converts standard ECG signals to VCG representations.  
   - Extracts critical ECG parameters.  
   - Computes electrical axes using  distinct definitions.  

3. **Polar Coordinate Conversion**  
   - Transforms 3D vector data into polar coordinates for further analysis.
