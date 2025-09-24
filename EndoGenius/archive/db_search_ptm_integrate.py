# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 10:33:13 2024

@author: lafields2
"""

# from database_search import raw_file_detail_extraction

# raw_file_formatted_path = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_input\fasta_formatted\2021_0817_CoG_1.txt"
# output_parent_directory = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.1.2\05"
# raw_file_detail_extraction(raw_file_formatted_path,output_parent_directory)


# predefined_db_path = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.1.2\unit_test_input\small_w_dileu.csv"
# output_parent_directory = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.1.2\06\2021_0817_CoG_1"
# choose_mzml_directory = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_input\fasta_formatted\2021_0817_CoG_1.ms2"
# raw_file_formatted_path = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_input\db_formatted\2021_0817_CoG_1.txt"
# precursor_error_cutoff = 20
# fragment_error_cutoff = 0.02
# min_mz = 200
# min_intensity = 1000
# standard_err_percent = 0.1
# amidation = 1
# oxidation_M_status = 1
# pyroglu_E_status = 1
# pyroglu_Q_status = 1
# sulfo_Y_status = 1
# dileu_status = 1
# max_modifications = 2
# confident_seq_cov = 70
# max_adjacent_swapped_AAs = 2
# min_motif_len = 3
# max_swapped_AA = 1
# max_adjacent_swapped_AA = 2
# sample_output_directory = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.1.2\06\2021_0817_CoG_1"

# from database_search import launch_db_search_pt1
# launch_db_search_pt1(predefined_db_path,output_parent_directory,choose_mzml_directory,raw_file_formatted_path,precursor_error_cutoff,
#                          fragment_error_cutoff,min_mz,min_intensity,standard_err_percent,amidation,oxidation_M_status,
#                          pyroglu_E_status,pyroglu_Q_status,sulfo_Y_status,dileu_status,max_modifications,sample_output_directory)



from PSM_assignment import PSM_assignment_execute

# db_search_parent_directory = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.1.2\06"
# motif_db_path = r"D:\Manuscripts\2024_EndoGeniusDIA\OG_SL_files\EndoGenius_DDA_search_results\input_files\NP_motif_db.csv"
# sample_output_directory = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.1.2\06\2021_0817_CoG_1"
# target_path = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.1.2\unit_test_input\small_w_dileu.csv"

standard_err_percent = 0.1
confident_seq_cov = 70.0
max_adjacent_swapped_AA = 2
min_motif_len = 3
fragment_error_threshold = 0.02
num_sub_AAs = 1
db_search_parent_directory = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.1.2\v07"
target_path = r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.1.2\v07\target_list.csv"
motif_path = r"D:\Manuscripts\2024_EndoGeniusDIA\OG_SL_files\EndoGenius_DDA_search_results\input_files\NP_motif_db.csv"

PSM_assignment_execute(standard_err_percent,
                       confident_seq_cov,
                       max_adjacent_swapped_AA,
                       min_motif_len,
                       fragment_error_threshold,
                       num_sub_AAs,
                       db_search_parent_directory,
                       target_path,
                       motif_path,
                       sample_output_directory)
