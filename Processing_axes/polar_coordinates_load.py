#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 17:57:21 2024

@author: mka23
"""


import pandas as pd

# Load precomputed results from CSV
polar_coordinates_healthy = pd.read_csv('polar_coordinates_results/polar_coordinates_healthy.csv')
polar_coordinates_healthy_inv = pd.read_csv('polar_coordinates_results/polar_coordinates_healthy_inv.csv')
polar_coordinates_all = pd.read_csv('polar_coordinates_results/polar_coordinates_all.csv')
polar_coordinates_all_inv = pd.read_csv('polar_coordinates_results/polar_coordinates_all_inv.csv')
