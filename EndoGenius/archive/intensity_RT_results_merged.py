# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 16:57:20 2024

@author: lafields2
"""

import pandas as pd
import csv

raw_file_details_path = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.0.5\iteration4\2021_0817_CoG_1_formatted.txt"
peptide_score_path = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.0.5\iteration4\2021_0817_CoG_1\IDed_peptide_scores.csv"

raw_converter = pd.read_csv(raw_file_details_path, sep=",",skiprows=[0], names= ['fragment_mz',
                                                                              'fragment_resolution',
                                                                              'fragment_z',
                                                                              'fragment_intensity',
                                                                              'precursor_mz',
                                                                              'ms2_scan',
                                                                              'precursor_z',
                                                                              'precursor_RT',
                                                                              'IonInjectTime',
                                                                              'ms1_scan',
                                                                              'precursor_intensity',
                                                                              'null'])

raw_converter_filtered = raw_converter.drop_duplicates(subset='ms2_scan')

dsd_summary_results = pd.read_csv(peptide_score_path)

dsd_summary_results2 = dsd_summary_results.merge(raw_converter_filtered, left_on='Scan', right_on='ms2_scan', how='left')
dsd_summary_results2 = dsd_summary_results2.drop_duplicates()

dsd_summary_results3 = dsd_summary_results2.drop(columns={'fragment_mz',
                                                 'fragment_resolution',
                                                 'fragment_z',
                                                 'fragment_intensity',
                                                 'precursor_mz',
                                                 'ms2_scan',
                                                 'precursor_z',
                                                 'null'})



file_path = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.0.5" + '\\merge_test.csv'
with open(file_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        dsd_summary_results2.to_csv(filec,index=False)
        
file_path = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.0.5" + '\\merge_test2.csv'
with open(file_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        dsd_summary_results3.to_csv(filec,index=False)
