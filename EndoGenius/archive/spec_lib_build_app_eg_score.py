# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 11:15:15 2023

@author: lawashburn
"""

import pandas as pd
import csv

results_directory = r"D:\Manuscripts\2024_EndoGeniusDIA\OG_SL_files\EndoGenius_DDA_search_results\20190701_CoG_Unlabeled_DDA_TR3\20190701_CoG_Unlabeled_DDA_TR3"
formatted_spectra_path = r"D:\Manuscripts\2024_EndoGeniusDIA\OG_SL_files\EndoGenius_DDA_search_results\20190701_CoG_Unlabeled_DDA_TR3\20190701_CoG_Unlabeled_DDA_TR3_formatted.txt"
output_directory = r"D:\Manuscripts\2024_EndoGeniusDIA\OG_SL_files\EndoGenius_DDA_search_results\intermediate_spectral_libraries"
fragment_error = 0.02

backslash_index1 = results_directory.rfind('\\')
backslash_index2 = results_directory.rfind('/')

if backslash_index1 > backslash_index2:
    backslash_index = backslash_index1
elif backslash_index2 >= backslash_index1:
    backslash_index = backslash_index2

base_file_path = results_directory[0:(backslash_index+1)]

tissue_type = results_directory.replace(base_file_path,'')

final_target_results_path = results_directory + '\\final_results_EG_score.csv'
final_target_results = pd.read_csv(final_target_results_path)

sequence_log = []
mz_log = []
z_log = []
rt_log = []
peaks_count_log = []
peak_list_log = []
sample_name_log = []
score_log = []
prec_int_log = []
hyper_log = []
per_seq_cov_log = []

formatted_spectra = pd.read_csv(formatted_spectra_path, sep=",",skiprows=[0], names= ['fragment_mz',
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

for row in final_target_results.index:
    peak_summary = []
    
    sequence = final_target_results['Peptide'][row]
    scan = final_target_results['Scan'][row]
    final_score = final_target_results['Final score, run: 1'][row]
    prec_int = final_target_results['precursor_intensity'][row]
    hyperscore = final_target_results['Hyperscore'][row]
    per_seq_cov = final_target_results['% Seq Coverage'][row]
    
    sequence_formatted = sequence.replace('(Glu->pyro-Glu)','(pyroGlu)')
    sequence_formatted = sequence_formatted.replace('(Gln->pyro-Glu)','(pyroGlu)')
    
    sample_fragment_path = results_directory + '\\fragment_matches\\' + sequence_formatted + '_' + str(scan) + '_fragment_report.ftr'
    
    sample_fragment = pd.read_feather(sample_fragment_path)
    
    sample_fragment_filtered = sample_fragment[sample_fragment['Fragment error (Da)'] <= fragment_error]

    spectra_scan_filtered = formatted_spectra[formatted_spectra['ms2_scan'] == scan]
    rt = spectra_scan_filtered['precursor_RT'].iloc[0]

    for line in sample_fragment_filtered.index:

        frag_name = sample_fragment_filtered['ion'][line]
        frag_mz = sample_fragment_filtered['Fragment actual m/z'][line]
        frag_z = sample_fragment_filtered['Fragment actual charge'][line]
        frag_int = sample_fragment_filtered['Fragment actual intensity'][line]        
        summary = frag_name + ':' + str(frag_mz) + ':' + str(frag_z) + ':' + str(frag_int)
        peak_summary.append(summary)

    
    mz = sample_fragment_filtered['Precursor actual m/z'].iloc[0]
    z = sample_fragment_filtered['Precursor actual charge'].iloc[0]
    length = len(sample_fragment_filtered)
    
    sequence_log.append(sequence)
    mz_log.append(mz)
    z_log.append(z)
    peaks_count_log.append(length)
    peak_list_log.append(peak_summary)
    rt_log.append(rt)
    score_log.append(final_score)
    sample_name_log.append(tissue_type)
    prec_int_log.append(prec_int)
    hyper_log.append(hyperscore)
    per_seq_cov_log.append(per_seq_cov)
    
individual_library = pd.DataFrame(
    {'Sequence': sequence_log,
     'Sample':sample_name_log,
     'EndoGenius score':score_log,
     'precursor m/z': mz_log,
     'precursor z': z_log,
     'precursor RT':rt_log,
     'precursor intensity':prec_int_log,
     'hyperscore':hyper_log,
     'fragment peak count': peaks_count_log,
     'fragment peak list':peak_list_log
    })

individual_library['Norm. hyperscore'] = individual_library['hyperscore'] / max(individual_library['hyperscore'])
individual_library['Norm. precursor intensity'] = individual_library['precursor intensity'] / max(individual_library['precursor intensity'])
individual_library['Norm. EndoGenius score'] = individual_library['EndoGenius score'] / max(individual_library['EndoGenius score'])
individual_library['Norm. precursor intensity'] = individual_library['precursor intensity'] / max(individual_library['precursor intensity'])

output_path = output_directory + '\\spectral_library_' + tissue_type + '.csv'
with open(output_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        individual_library.to_csv(filec,index=False)
