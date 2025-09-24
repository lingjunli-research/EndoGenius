# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 16:10:11 2024

@author: lafields2
"""

import pandas as pd
import csv

quant_output_directory = r"D:\Manuscripts\2024_MotifAutomation\Algorithm_Validation\EndoGenius_search_results\NP_search_of_AMP_hemolymph"

quant_file_paths = [r"D:\Manuscripts\2024_MotifAutomation\Algorithm_Validation\EndoGenius_search_results\NP_search_of_AMP_hemolymph\BR1_TR1\Control1_BR1_TR1\final_results_EG_score.csv",
                    r"D:\Manuscripts\2024_MotifAutomation\Algorithm_Validation\EndoGenius_search_results\NP_search_of_AMP_hemolymph\BR1_TR2\Control1_BR1_TR2\final_results_EG_score.csv",
                    r"D:\Manuscripts\2024_MotifAutomation\Algorithm_Validation\EndoGenius_search_results\NP_search_of_AMP_hemolymph\BR1_TR3\Control1_BR1_TR3\final_results_EG_score.csv",
                    r"D:\Manuscripts\2024_MotifAutomation\Algorithm_Validation\EndoGenius_search_results\NP_search_of_AMP_hemolymph\BR2_TR1\Control2_BR2_TR1\final_results_EG_score.csv",
                    r"D:\Manuscripts\2024_MotifAutomation\Algorithm_Validation\EndoGenius_search_results\NP_search_of_AMP_hemolymph\BR2_TR2\Control2_BR2_TR2\final_results_EG_score.csv",
                    r"D:\Manuscripts\2024_MotifAutomation\Algorithm_Validation\EndoGenius_search_results\NP_search_of_AMP_hemolymph\BR2_TR3\Control2_BR2_TR3\final_results_EG_score.csv",
                    r"D:\Manuscripts\2024_MotifAutomation\Algorithm_Validation\EndoGenius_search_results\NP_search_of_AMP_hemolymph\BR3_TR1\Control3_BR3_TR1\final_results_EG_score.csv",
                    r"D:\Manuscripts\2024_MotifAutomation\Algorithm_Validation\EndoGenius_search_results\NP_search_of_AMP_hemolymph\BR3_TR2\Control3_BR3_TR2\final_results_EG_score.csv",
                    r"D:\Manuscripts\2024_MotifAutomation\Algorithm_Validation\EndoGenius_search_results\NP_search_of_AMP_hemolymph\BR3_TR3\Control3_BR3_TR3\final_results_EG_score.csv"]

quant_file_names = ['BR1_TR1',
                    'BR1_TR2',
                    'BR1_TR3',
                    'BR2_TR1',
                    'BR2_TR2',
                    'BR2_TR3',
                    'BR3_TR1',
                    'BR3_TR2',
                    'BR3_TR3']

quant_merge_df = pd.DataFrame()

if len(quant_file_paths) != len(quant_file_names):
    raise Exception("A file name must be selected for each file input")
    
else:
    for a in range(0,(len(quant_file_paths))):
        quant_single_file_path = quant_file_paths[a]
        quant_single_name = quant_file_names[a]
        
        quant_single_file = pd.read_csv(quant_single_file_path)
        quant_single_file_filtered = pd.DataFrame()
        quant_single_file_filtered['Peptide'] = quant_single_file['Peptide']
        quant_single_file_filtered[(quant_single_name + ' Intensity')] = quant_single_file['precursor_intensity']
        
        quant_single_file_filtered = quant_single_file_filtered.sort_values(by=(quant_single_name + ' Intensity'),ascending=False)
        quant_single_file_filtered = quant_single_file_filtered.drop_duplicates(subset='Peptide')

        if len(quant_merge_df)==0:
            quant_merge_df = quant_single_file_filtered
        else:
            quant_merge_df = pd.merge(quant_merge_df, quant_single_file_filtered, on='Peptide', how='outer')

quant_file_out_path = quant_output_directory + '\\merged_intensities.csv'
with open(quant_file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        quant_merge_df.to_csv(filec,index=False)
