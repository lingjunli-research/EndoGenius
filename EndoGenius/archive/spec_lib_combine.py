# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 16:55:21 2024

@author: lafields2
"""

import pandas as pd
import csv

list_of_files = [r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.0.8\sl_build_output\v02\spectral_library_2021_0817_CoG_1.csv",
                 r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.0.8\sl_build_output\v02\spectral_library_2021_0817_CoG_2.csv",
                 r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.0.8\sl_build_output\v02\spectral_library_2021_0817_CoG_3.csv"]

output_directory = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.0.8\sl_build_output\v02"

sl_w_duplicates = pd.DataFrame()

for file in list_of_files:
    spectra = pd.read_csv(file)
    sl_w_duplicates = pd.concat([sl_w_duplicates, spectra], ignore_index=True, sort=False)
    
output_path = output_directory + '\\combined_library_unfiltered.csv'
with open(output_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        sl_w_duplicates.to_csv(filec,index=False)

sl_filtered = sl_w_duplicates.sort_values(by='Score',ascending=False)
sl_filtered = sl_filtered.drop_duplicates(subset='Sequence')

sequence_storage = []
sample_storage = []

for row in sl_filtered.index:
    sequence = sl_filtered['Sequence'][row]
    
    sl_dups_filtered = sl_w_duplicates[sl_w_duplicates['Sequence'] == sequence]
    sl_dups_filtered = sl_dups_filtered.drop_duplicates(subset='Sample')
    
    sample_list = sl_dups_filtered['Sample'].tolist()
    
    sequence_storage.append(sequence)
    sample_storage.append(sample_list)

sequence_sample_log = pd.DataFrame(
    {'Sequence': sequence_storage,
     'Sample':sample_storage})

sl_filtered = sl_filtered.drop(columns=['Sample'])

sl_filtered = pd.merge(sl_filtered, sequence_sample_log, on='Sequence', how='outer')

output_path = output_directory + '\\combined_library_filtered.csv'
with open(output_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        sl_filtered.to_csv(filec,index=False)
