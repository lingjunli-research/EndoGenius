# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 17:56:16 2024

@author: lafields2
"""

import pandas as pd
import csv

results_path = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.1.2\v10\20240523_DiLeu_TR1_240523184202_5mod_5FDR\20240523_DiLeu_TR1_240523184202\final_results__target.csv"
spectra_path = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.1.2\v10\20240523_DiLeu_TR1_240523184202_5mod_5FDR\20240523_DiLeu_TR1_240523184202_formatted.txt"

output_path = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.1.4"

mass_h = 1.00784
fragment_error_threshold = 20

spectra = pd.read_csv(spectra_path, sep=",",skiprows=[0], names= ['fragment_mz',
                                                                              'fragment_intensity',
                                                                              'fragment_z',
                                                                              'fragment_resolution',
                                                                              'precursor_mz',
                                                                              'ms2_scan',
                                                                              'precursor_z',
                                                                              'precursor_RT',
                                                                              'IonInjectTime',
                                                                              'ms1_scan',
                                                                              'precursor_intensity',
                                                                              'null'])

spectra['fragment_z'] = spectra['fragment_z'].replace(0, 1)
spectra['fragment_monoisotopic_mass'] = (spectra['fragment_mz'] * spectra['fragment_z']) - (mass_h * spectra['fragment_z'])

results = pd.read_csv(results_path)

tag_name_list = ['115a','115b','116a','116b','116c','117a','117b','117c','118a','118b','118c','118d']
tag_mass_list = [115.12476,115.13108,
                 116.12812,116.13444,116.14028,
                 117.13147,117.13731,117.14363,
                 118.13483,118.14067,118.14699,118.15283]

spectra_w_details = pd.DataFrame()

peptide_name_log = []
tag_name_log = []
tag_intensity_log = []
scan_log = []
fragment_mz_log = []

for y in tag_mass_list:
    y_monoisotopic = (y*1) - (mass_h*1)
    spectra[str(y) + ' error'] = ((abs(spectra['fragment_monoisotopic_mass'] - y_monoisotopic))/y_monoisotopic)*1E6
    spectra_w_details = spectra

for a in range(0,len(results)):
    peptide_id = results['Peptide'].iloc[a]
    if '(12PlexDiLeu)' in peptide_id:
        scan_id = results['Scan'].iloc[a]

        spectra_filtered = spectra_w_details[spectra_w_details['ms2_scan'] == scan_id]
        
        for k in range(0,len(tag_name_list)):
            tag_name_selected = tag_name_list[k]
            tag_mass_selected = tag_mass_list[k]
            
            if spectra_filtered[str(tag_mass_selected) + ' error'].min() <= fragment_error_threshold:
                spectra_filtered_tag = spectra_filtered[spectra_filtered[str(tag_mass_selected) + ' error'] <= fragment_error_threshold]
                tag_intensity_log.append(spectra_filtered_tag['fragment_intensity'].max())
                tag_name_log.append(tag_name_selected)
                peptide_name_log.append(peptide_id)
                scan_log.append(scan_id)
                
                filter_filter_df = spectra_filtered_tag[spectra_filtered_tag['fragment_intensity'] == (spectra_filtered_tag['fragment_intensity'].max())]
                fragment_mz_log.append(filter_filter_df['fragment_mz'].iloc[0])
                
            else:
                tag_intensity_log.append(0)
                tag_name_log.append(tag_name_selected)
                peptide_name_log.append(peptide_id)
                scan_log.append(scan_id)
                fragment_mz_log.append(0)
        
report_ion_intensity_df = pd.DataFrame()
report_ion_intensity_df['Peptide'] = peptide_name_log
report_ion_intensity_df['Scan'] = scan_log
report_ion_intensity_df['Tag'] = tag_name_log
report_ion_intensity_df['Intensity'] = tag_intensity_log
report_ion_intensity_df['Fragment_mz'] = fragment_mz_log

# Pivot the table without hardcoding the values in the 'Tag' column
pivot_df = report_ion_intensity_df.pivot_table(index=['Peptide', 'Scan'], columns='Tag', values='Intensity').reset_index()

# Flatten the MultiIndex columns
pivot_df.columns.name = None
pivot_df.columns = ['Peptide', 'Scan'] + [f'Intensity_{tag}' for tag in pivot_df.columns[2:]]


output_path_rep = output_path + '\\reporter_ions_extracted.csv'
# df.reset_index().to_feather(output_path_rep) 

with open(output_path_rep,'w',newline='') as filec:
        writerc = csv.writer(filec)
        pivot_df.to_csv(filec,index=False)       
        
        
        
