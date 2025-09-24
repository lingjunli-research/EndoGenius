# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 12:27:33 2023

@author: lawashburn
"""

import pandas as pd
import numpy as np
import csv
import smtplib

print('Target-Decoy Assess')

def target_decoy_apply(dsd_summary_results,target_results,output_directory,fdr_cutoff,sample_output_directory,raw_file_formatted_path):

    number_runs = 1
    check_iterations = 5
    
    raw_converter = pd.read_csv(raw_file_formatted_path, sep=",",skiprows=[0], names= ['fragment_mz',
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
    target_IDS = pd.read_csv(target_results)
    dsd_summary_results = dsd_summary_results

    
    target_list_str = target_IDS['Sequence'].values.tolist()
    
    dsd_summary_results['Unmodified sequence'] = dsd_summary_results['Peptide']
    dsd_summary_results['Unmodified sequence'] = dsd_summary_results['Unmodified sequence'].str.replace(r'\s*\(.*?\)\s*', '', regex=True)
    dsd_summary_results['Status'] = dsd_summary_results['Unmodified sequence'].apply(lambda x: any([k in x for k in target_list_str]))

    samples_df = dsd_summary_results.drop_duplicates(subset='Sample')
    
    
    
    samples = samples_df['Sample'].values.tolist()
    
    for sample in samples:
        
        run_log = []
        num_target_IDs_log = []
        num_unique_IDs = []
        num_decoy_IDs = []
        
        dsd_summary_results_per_sample_pre = dsd_summary_results[dsd_summary_results['Sample'] == sample]
        
        dsd_summary_results2 = dsd_summary_results_per_sample_pre.merge(raw_converter_filtered, left_on='Scan', right_on='ms2_scan', how='left')
        dsd_summary_results2 = dsd_summary_results2.drop_duplicates()
        
        dsd_summary_results_per_sample = dsd_summary_results2.drop(columns={'fragment_mz',
                                                         'fragment_resolution',
                                                         'fragment_z',
                                                         'fragment_intensity',
                                                         'precursor_mz',
                                                         'ms2_scan',
                                                         'precursor_z',
                                                         'null'})
        
        target_results = dsd_summary_results_per_sample[dsd_summary_results_per_sample['Status'] == True] 
        decoy_results = dsd_summary_results_per_sample[dsd_summary_results_per_sample['Status'] == False]
        for a in range(1,(number_runs+1)):
                target_scores = target_results['Final score, run: ' + str(a)].values.tolist()
                decoy_scores = decoy_results['Final score, run: ' + str(a)].values.tolist()
                
                all_scores = dsd_summary_results['Final score, run: ' + str(a)].values.tolist()
                min_scores = 0
                max_scores = (max(all_scores)) +10
                int_list_no_dups = list(set(target_scores))
                score_list_sorted = list(reversed(sorted(int_list_no_dups)))
                
                score_cutoff = []
                fdr = []
                num_target_IDs = []
                
                for score in score_list_sorted:
                    if len(fdr)>0:
                        if fdr[-1] <= fdr_cutoff*check_iterations:
                            target_results_filtered = target_results[target_results['Final score, run: ' + str(a)] >= score]
                            decoy_results_filtered = decoy_results[decoy_results['Final score, run: ' + str(a)] >= score]
                            if len(target_results_filtered) >0:
                                fdr_prelim = len(decoy_results_filtered)/len(target_results_filtered)
                                fdr_report = round(fdr_prelim,2)
                            else:
                                fdr_report = 1
                            
                            score_cutoff.append(score)
                            fdr.append(fdr_report)
                            if len(target_results_filtered)>0:
                                num_target_IDs.append(len(target_results_filtered))
                            else:
                                num_target_IDs.append(0)
                           
                            if len(decoy_results_filtered)>0:
                                num_decoy_IDs.append(len(decoy_results_filtered))
                            else:
                                num_decoy_IDs.append(0)
                        else:
                            pass
                    elif len(fdr) == 0:
                        target_results_filtered = target_results[target_results['Final score, run: ' + str(a)] >= score]
                        decoy_results_filtered = decoy_results[decoy_results['Final score, run: ' + str(a)] >= score]
                        if len(target_results_filtered) >0:
                            fdr_report = len(decoy_results_filtered)/len(target_results_filtered)
                        else:
                            fdr_report = 1
                        
                        score_cutoff.append(score)
                        fdr.append(fdr_report)
                        if len(target_results_filtered)>0:
                            num_target_IDs.append(len(target_results_filtered))
                        else:
                            num_target_IDs.append(0)
                       
                        if len(decoy_results_filtered)>0:
                            num_decoy_IDs.append(len(decoy_results_filtered))
                        else:
                            num_decoy_IDs.append(0)
                    else:
                        pass
                
                fdr_table = pd.DataFrame()
                fdr_table['Score Threshold'] = score_cutoff
                fdr_table['FDR'] = fdr
                fdr_table['# Target IDs'] = num_target_IDs
                
                fdr_table_format = fdr_table
                fdr_table_format['FDR'] = fdr_table_format['FDR']*100
                fdr_table_format = fdr_table_format.rename(columns={'FDR':'FDR (%)','Score Threshold':'EndoGenius Score'})
                fdr_table_format['log(EndoGenius Score)'] = np.log10(fdr_table_format['EndoGenius Score'])
                fig = fdr_table_format.plot.scatter(x = 'FDR (%)', y = 'log(EndoGenius Score)',c='blue').get_figure()
                fig_out_path = sample_output_directory + '\\fdr_v_score.svg'
                fig.savefig(fig_out_path, dpi=1500)
                
                file_path = sample_output_directory +'\\FDR_eval_table.csv'
                with open(file_path,'w',newline='') as filec:
                        writerc = csv.writer(filec)
                        fdr_table.to_csv(filec,index=False)
                        
                fdr_results_filtered = fdr_table[fdr_table['FDR'] <= fdr_cutoff]
                
                best_num_targets = fdr_results_filtered['# Target IDs'].max()
                
                unique_peptide_score_thresh_filter = fdr_table[fdr_table['# Target IDs'] == best_num_targets]
                unique_peptide_score_thresh = unique_peptide_score_thresh_filter['Score Threshold'].min()
                
                filter_target_score_thresh = target_results[target_results['Final score, run: ' + str(a)] >= unique_peptide_score_thresh]
                filter_target_score_thresh_unique = filter_target_score_thresh.drop_duplicates(subset=['Peptide'])
                
                run_log.append(a)
                num_target_IDs_log.append(best_num_targets)
                num_unique_IDs.append(len(filter_target_score_thresh_unique))
                
                
                fdr_table_filtered = fdr_table[fdr_table['# Target IDs'] == best_num_targets]
                score_threshold_to_apply = fdr_table_filtered['Score Threshold'].min()
                
                target_results_final = target_results[target_results['Final score, run: ' + str(a)] >= score_threshold_to_apply]
                decoy_results_final = decoy_results[decoy_results['Final score, run: ' + str(a)] >= score_threshold_to_apply]
                
                file_path = sample_output_directory + '\\final_results__decoy.csv'
                with open(file_path,'w',newline='') as filec:
                        writerc = csv.writer(filec)
                        decoy_results_final.to_csv(filec,index=False)
                
                file_path = sample_output_directory + '\\final_results__target.csv'
                with open(file_path,'w',newline='') as filec:
                        writerc = csv.writer(filec)
                        target_results_final.to_csv(filec,index=False)
            
        dsd_merge_table = pd.DataFrame()
        dsd_merge_table['Run #'] = run_log
        dsd_merge_table['# Target IDs'] = num_target_IDs_log
        dsd_merge_table['# Unique Target IDs'] = num_unique_IDs
        return dsd_merge_table