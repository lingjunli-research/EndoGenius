# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 12:27:33 2023

@author: lawashburn
"""

import pandas as pd
import numpy as np
import csv
import smtplib

def target_decoy_apply(dsd_summary_results,target_results,output_directory,fdr_cutoff):
    #dsd_summary_results = r"C:\Users\lawashburn\Documents\DB_pep_validation\motif_percent_assess_2023010\DSD_results_scores.csv"
    #target_results = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\Reference_DB\target_list.csv"
    #output_directory = r"C:\Users\lawashburn\Documents\DB_pep_validation\motif_percent_assess_2023010\100_percent_bulk"
    
    #fdr_cutoff = 0.01
    number_runs = 1
    check_iterations = 5000
    
    
    target_IDS = pd.read_csv(target_results)
    dsd_summary_results = dsd_summary_results
    
    target_list_str = target_IDS['Sequence'].values.tolist()
    
    dsd_summary_results['Unmodified sequence'] = dsd_summary_results['Peptide']
    dsd_summary_results['Unmodified sequence'] = dsd_summary_results['Unmodified sequence'].str.replace(r"\([^)]*\)","")
    dsd_summary_results['Status'] = dsd_summary_results['Unmodified sequence'].apply(lambda x: any([k in x for k in target_list_str]))
    
    
    
    
    
    samples_df = dsd_summary_results.drop_duplicates(subset='Sample')
    samples = samples_df['Sample'].values.tolist()
    
    for sample in samples:
        
        run_log = []
        num_target_IDs_log = []
        num_unique_IDs = []
        num_decoy_IDs = []
        
        dsd_summary_results_per_sample = dsd_summary_results[dsd_summary_results['Sample'] == sample]
        target_results = dsd_summary_results_per_sample[dsd_summary_results_per_sample['Status'] == True]
        decoy_results = dsd_summary_results_per_sample[dsd_summary_results_per_sample['Status'] == False]
        #for ind in dsd_summary_results_per_sample.index:
        for a in range(1,(number_runs+1)):
                target_scores = target_results['Final score, run: ' + str(a)].values.tolist()
                decoy_scores = decoy_results['Final score, run: ' + str(a)].values.tolist()
                
                all_scores = dsd_summary_results['Final score, run: ' + str(a)].values.tolist()
                min_scores = 0
                max_scores = (max(all_scores)) +10
                
                #bins = np.linspace(min_scores, max_scores, 50)
                        
                int_list_no_dups = list(set(target_scores))
                
                # for aa in target_scores:
                #     aaa = int(aa)
                #     if aaa not in int_list_no_dups:
                #         int_list_no_dups.append(aaa)
                
                
                score_list_sorted = list(reversed(sorted(int_list_no_dups)))
                
                score_cutoff = []
                fdr = []
                num_target_IDs = []
                
                for score in score_list_sorted:
                    if len(fdr)>0:
                        if fdr[-1] <= fdr_cutoff*5:
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
                
                file_path = output_directory + '\\FDR_eval_table_' +  sample + '_DSD_run_'+ str(a) + '.csv'
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
                
                file_path = output_directory + '\\final_results_' +  sample + '_DSD_run_'+ str(a) + '_decoy.csv'
                with open(file_path,'w',newline='') as filec:
                        writerc = csv.writer(filec)
                        decoy_results_final.to_csv(filec,index=False)
                
                file_path = output_directory + '\\final_results_' +  sample + '_DSD_run_'+ str(a) + '_target.csv'
                with open(file_path,'w',newline='') as filec:
                        writerc = csv.writer(filec)
                        target_results_final.to_csv(filec,index=False)
            
        dsd_merge_table = pd.DataFrame()
        dsd_merge_table['Run #'] = run_log
        dsd_merge_table['# Target IDs'] = num_target_IDs_log
        dsd_merge_table['# Unique Target IDs'] = num_unique_IDs
        
            
        file_path = output_directory + '\\dsd_eval_' +  sample + '.csv'
        with open(file_path,'w',newline='') as filec:
                writerc = csv.writer(filec)
                dsd_merge_table.to_csv(filec,index=False)
        
        return dsd_merge_table
                    
        # subject = 'Your code has finished running'
        # text = 'FDR-DSD evaluation completed for ' + sample
        # content = 'Subject: %s\n\n%s' % (subject, text)
        # mail = smtplib.SMTP('smtp.gmail.com',587)
        # mail.ehlo()
        # mail.starttls()
        # mail.login('lingjun.li.notifications@gmail.com','eabtnjwaikdssdtd')
        # mail.sendmail('lingjun.li.notifications@gmail.com','lawashburn@wisc.edu',content) 
        # mail.close()
    
    subject = 'Your code has finished running'
    text = 'FDR-DSD evaluation completed for all samples'
    content = 'Subject: %s\n\n%s' % (subject, text)
    mail = smtplib.SMTP('smtp.gmail.com',587)
    mail.ehlo()
    mail.starttls()
    mail.login('lingjun.li.notifications@gmail.com','eabtnjwaikdssdtd')
    mail.sendmail('lingjun.li.notifications@gmail.com','lawashburn@wisc.edu',content) 
    mail.close()
