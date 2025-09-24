# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 09:11:25 2024

@author: lafields2
"""

import pandas as pd
import csv

output_directory = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.0.7\iteration2\output"

file_paths = [r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.0.7\iteration2\input\final_results_target_CoG1.csv",
              r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.0.7\iteration2\input\final_results_target_CoG2.csv",
              r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.0.7\iteration2\input\final_results_target_CoG3.csv"]

file_names = ['CoG1','CoG2','CoG3']

merge_df = pd.DataFrame()

if len(file_paths) != len(file_names):
    raise Exception("A file name must be selected for each file input")
    
else:
    for a in range(0,(len(file_paths))):
        single_file_path = file_paths[a]
        single_name = file_names[a]
        
        single_file = pd.read_csv(single_file_path)
        single_file_filtered = pd.DataFrame()
        single_file_filtered['Peptide'] = single_file['Peptide']
        single_file_filtered[(single_name + ' Intensity')] = single_file['precursor_intensity']
        
        single_file_filtered = single_file_filtered.sort_values(by=(single_name + ' Intensity'),ascending=False)
        single_file_filtered = single_file_filtered.drop_duplicates(subset='Peptide')

        if len(merge_df)==0:
            merge_df = single_file_filtered
        else:
            merge_df = pd.merge(merge_df, single_file_filtered, on='Peptide', how='outer')

file_out_path = output_directory + '\\merged_intensities.csv'
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        merge_df.to_csv(filec,index=False)
