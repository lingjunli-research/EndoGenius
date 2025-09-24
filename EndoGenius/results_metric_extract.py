# -*- coding: utf-8 -*-
"""
Created on Tue May  9 14:30:59 2023

@author: lawashburn
"""
import os
import pandas as pd
import re
import csv
from statistics import mean
import smtplib

print('Results Metric Extract')

def results_metric_extract_start(results_directory,output_directory,sample_output_directory):

    def get_dir_names_with_strings_list(str_list): #definition for finding a file containing a string in filename in specified directory
        full_list = [ name for name in os.listdir(results_directory) if os.path.isdir(os.path.join(results_directory, name)) ]
        final_list = [nm for ps in str_list for nm in full_list if ps in nm]
        return final_list
    
    def get_file_names_with_strings_list(str_list, folder): #definition for finding a file containing a string in filename in specified directory
        full_list = os.listdir(folder)
        final_list = [nm for ps in str_list for nm in full_list if ps in nm]
        return final_list
    
    def count_consec(lst):
        consec = [1]
        for x, y in zip(lst, lst[1:]):
            if x == y - 1:
                consec[-1] += 1
            else:
                consec.append(1)
        return consec
    
    query = '' #search for ion list pertaining to the sequence
    parent_dir_list = (get_dir_names_with_strings_list([query])) #search for the file based on query
    
    peptide_storage = []
    scan_storage = []
    sample_storage = []
    avg_fragment_err_storage = []
    precursor_err_storage = []
    motif_score_storage = []
    consec_b_ions_storage = []
    consec_y_ions_storage = []
    c_term_amid_storage = []
    percent_seq_cov_storage = []
    avg_annotate_peak_storage = []
    percent_b_ions_storage = []
    percent_y_ions_storage = []
    avg_ions_per_AA_storage = []
    num_non_neut_ions_storage = []
    hyperscore_storage = []
    e_pyroglu = []
    q_pyroglu = []
    max_fragment_err_storage = []
    min_fragment_err_storage = []
    med_fragment_err_storage = []
    ox_storage = []
    sulfo_storage = []
    num_mods_storage = []
    seq_cov_storage = []

    for a in parent_dir_list:
        a_path = sample_output_directory
        csv_list = (get_file_names_with_strings_list(['.csv'],a_path)) #search for the file based on query
    
        final_results_file_path = sample_output_directory + '\\psm_results_w_motif.csv'
        final_results = pd.read_csv(final_results_file_path)
        final_results = final_results[final_results['Sequence coverage'] > 0]
        
        for ind in final_results.index:
            sequence = final_results['Sequence with mod'][ind]
            scan = int(final_results['Scan'][ind])
            motif_score = final_results['Final motif score'][ind]
            corr_score = final_results['Correlation value'][ind]
            seq_cov = final_results['Sequence coverage'][ind]
            
            peptide_storage.append(sequence)
            scan_storage.append(scan)
            motif_score_storage.append(motif_score)
            sample_storage.append(a)
            hyperscore_storage.append(corr_score)
            seq_cov_storage.append(seq_cov)
            
            num_mods = sequence.count('(')
            num_mods_storage.append(num_mods)
            
            if '(Amidated)' in sequence:
                c_term_amid_storage.append(1)
            else:
                c_term_amid_storage.append(0)
                
            if '(Oxidation)' in sequence:
                ox_storage.append(1)
            else:
                ox_storage.append(0)
                
            if '(Sulfo)' in sequence:
                sulfo_storage.append(1)
            else:
                sulfo_storage.append(0)
            
            if 'Q(Gln->pyro-Glu)' in sequence:
                q_pyroglu.append(1)
                e_pyroglu.append(0)
                updated_seq = sequence.replace('Q(Gln->pyro-Glu)','Q(pyroGlu)')
                fragment_report_path = sample_output_directory + '\\fragment_matches\\' + updated_seq + '_' + str(scan) + '_fragment_report.csv'
                
            elif 'E(Glu->pyro-Glu)' in sequence:
                e_pyroglu.append(1)
                q_pyroglu.append(0)
                updated_seq = sequence.replace('E(Glu->pyro-Glu)','E(pyroGlu)')
                fragment_report_path = sample_output_directory + '\\fragment_matches\\' + updated_seq + '_' + str(scan) + '_fragment_report.csv'
            
            else:
                fragment_report_path = sample_output_directory + '\\fragment_matches\\' + sequence + '_' + str(scan) + '_fragment_report.csv'
                e_pyroglu.append(0)
                q_pyroglu.append(0)
            
            fragment_report = pd.read_csv(fragment_report_path)
            fragment_report['Fragment error (Da)'] = fragment_report['Fragment error (Da)'].abs()
            fragment_err_mean = fragment_report['Fragment error (Da)'].mean()
            avg_fragment_err_storage.append(fragment_err_mean)
            fragment_err_max = fragment_report['Fragment error (Da)'].max()
            max_fragment_err_storage.append(fragment_err_max)
            fragment_err_min = fragment_report['Fragment error (Da)'].min()
            min_fragment_err_storage.append(fragment_err_min)
            fragment_err_med = fragment_report['Fragment error (Da)'].median()
            med_fragment_err_storage.append(fragment_err_med)
            
            precursor_err = fragment_report['Precursor error (ppm)_x'][0]
            precursor_err_storage.append(precursor_err)
            
            fragment_ion_list = fragment_report['ion'].values.tolist()
            num_ions = len(fragment_ion_list)
            
            b_ions = []
            y_ions = []
            
            for ion in fragment_ion_list:
                if 'b' in ion:
                    b_ions.append(ion)
                elif 'y' in ion:
                    y_ions.append(ion)
            
            percent_b_ions = (len(b_ions) / num_ions)*100
            percent_y_ions = (len(y_ions) / num_ions)*100
            
            percent_b_ions_storage.append(percent_b_ions)
            percent_y_ions_storage.append(percent_y_ions)
    
            peptide_length = len(sequence)
            ion_types = []       
            for i in range(2,peptide_length+1):
                ion_types.append('b'+str(i))        
            for i in range(1,peptide_length):
                ion_types.append('y'+str(i))            
            ion_list_present = fragment_report['ion'].values.tolist()
            ion_format_list = []        
            for aa in ion_list_present:
                aa = aa.replace('-H2O','')
                aa = aa.replace('-NH3','')
                ion_format_list.append(aa)        
            ion_counts = []        
            for ion in ion_types:
                count = ion_format_list.count(ion)
                ion_counts.append(count)        
            mean_ion_per_peak = mean(ion_counts)
            avg_annotate_peak_storage.append(mean_ion_per_peak)
    
            b_ion_loc = []
            for ion in b_ions:
                ion = ion.replace('b','')
                ion = ion.replace('-H2O','')
                ion = ion.replace('-NH3','')
                b_ion_loc.append(int(ion))
            
            y_ion_loc = []
            for ion in y_ions:
                ion = ion.replace('y','')
                ion = ion.replace('-H2O','')
                ion = ion.replace('-NH3','')
                y_ion_loc.append(int(ion))
            
            peptide_no_mods = re.sub("[\(\[].*?[\)\]]", "", sequence)
            num_expected_ions = len(peptide_no_mods)
            
            consec_b_ion = count_consec(b_ion_loc)
            max_consec_b_ion = max(consec_b_ion)
            norm_consec_b_ion = max_consec_b_ion / num_expected_ions
            consec_b_ions_storage.append(norm_consec_b_ion)
            
            consec_y_ion = count_consec(y_ion_loc)
            max_consec_y_ion = max(consec_y_ion)
            norm_consec_y_ion = max_consec_y_ion / num_expected_ions
            consec_y_ions_storage.append((norm_consec_y_ion))
            
            
            
            y_reorient = []
            
            for ion_loc in y_ion_loc:
                new_loc = num_expected_ions-ion_loc
                y_reorient.append(new_loc)
            
            all_ions = []
            all_ions_no_dups = []
            
            for i in b_ion_loc:
                all_ions.append(i)
                if i not in all_ions_no_dups:
                    all_ions_no_dups.append(i)
            
            for i in y_reorient:
                all_ions.append(i)
                if i not in all_ions_no_dups:
                    all_ions_no_dups.append(i)
            
            counts_per_aa = []
            
            for k in all_ions_no_dups:
                ammount = all_ions.count(k)
                counts_per_aa.append(ammount)
            
            mean_counts = mean(counts_per_aa)
            avg_ions_per_AA_storage.append(mean_counts)
            
    
            ion_filtered_loc = []
            ion_filtered_loc_no_dups = []
            
            for ion_nonfiltered in fragment_ion_list:
                if '-H2O' in ion_nonfiltered:
                    pass
                elif '-NH3' in ion_nonfiltered:
                    pass
                else:
                    if 'b' in ion_nonfiltered:
                        ion_formatted = ion_nonfiltered.replace('b','')
                        ion_filtered_loc.append(ion_formatted)
                        if ion_formatted not in ion_filtered_loc_no_dups:
                            ion_filtered_loc_no_dups.append(ion_formatted)
                        else:
                            pass
                    elif 'y' in ion_nonfiltered:
                        ion_formatted = ion_nonfiltered.replace('y','')
                        ion_filtered_loc.append(ion_formatted)
                        if ion_formatted not in ion_filtered_loc_no_dups:
                            ion_filtered_loc_no_dups.append(ion_formatted)
                        else:
                            pass
                    else:
                        pass
            
            filtered_counts_per_aa = []
            
            for k in ion_filtered_loc_no_dups:
                ammount = ion_filtered_loc.count(k)
                filtered_counts_per_aa.append(ammount)
            
            if len(filtered_counts_per_aa)>0:
                filtered_mean_counts = mean(filtered_counts_per_aa)
            else:
                filtered_mean_counts = 0
            num_non_neut_ions_storage.append(filtered_mean_counts)
    
    final_weighting_metrics = pd.DataFrame()
    final_weighting_metrics['Sample'] = sample_storage
    final_weighting_metrics['Peptide'] = peptide_storage
    final_weighting_metrics['Scan'] = scan_storage
    final_weighting_metrics['Average Fragment Error'] = avg_fragment_err_storage
    final_weighting_metrics['Precursor Error'] = precursor_err_storage 
    final_weighting_metrics['# Consecutive b-ions'] = consec_b_ions_storage
    final_weighting_metrics['# Consecutive y-ions'] = consec_y_ions_storage
    final_weighting_metrics['C-termini amidation'] = c_term_amid_storage
    final_weighting_metrics['Motif_Score'] = motif_score_storage
    final_weighting_metrics['Average annotations/fragment '] = avg_annotate_peak_storage
    final_weighting_metrics['% Fragment ions are b'] = percent_b_ions_storage
    final_weighting_metrics['% Fragment ions are y'] = percent_y_ions_storage
    final_weighting_metrics['Average number of fragment ions per AA'] = avg_ions_per_AA_storage
    final_weighting_metrics['Number of non-neutral-loss fragment ions per AA'] = num_non_neut_ions_storage
    final_weighting_metrics['Hyperscore'] = hyperscore_storage
    final_weighting_metrics['Pyro-glu on E'] = e_pyroglu
    final_weighting_metrics['Pyro-glu on Q'] = q_pyroglu
    final_weighting_metrics['Max Fragment Error'] = max_fragment_err_storage
    final_weighting_metrics['Min Fragment Error'] = min_fragment_err_storage
    final_weighting_metrics['Median Fragment Error'] = med_fragment_err_storage
    final_weighting_metrics['Oxidation'] = ox_storage
    final_weighting_metrics['Sulfation'] = sulfo_storage
    final_weighting_metrics['# Modifications'] = num_mods_storage
    final_weighting_metrics['% Seq Coverage'] = seq_cov_storage
    
    file_path = sample_output_directory + '\\metrics_extracted_motif.csv'
    with open(file_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            final_weighting_metrics.to_csv(filec,index=False)
    
    return final_weighting_metrics

