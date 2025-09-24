# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 10:02:33 2023

@author: lawashburn
"""

import pandas as pd
import numpy as np

import csv

print('Metric Handling')

def metric_handling_apply(metrics_extracted_path,output_directory,sample_output_directory):

    metrics = metrics_extracted_path
    metrics = metrics.drop(columns=['Pyro-glu on E','Pyro-glu on Q','Oxidation','Sulfation', '# Modifications','C-termini amidation',
                                    'Max Fragment Error','Min Fragment Error','Median Fragment Error','% Fragment ions are b',
                                    '% Fragment ions are y'])
    
    metrics_w_dsd_results = metrics
        
    run = 1
    precursor_err_mult = 10
    no_consec_b_ions_mult = 10
    seq_cov_mult = 10
    non_neut_frag_ions_per_AA_mult = 10
    consec_b_fact = 'Multiply'
    percent_seq_cov_fact = 'Multiply'
    non_neut_per_AA_fact = 'Divide'
    prec_err_fact = 'Divide'
    
    final_score_storage = []
    
    for b in range(0,len(metrics)):
    
        metrics_results = metrics.iloc[[b]]
        sample = metrics_results['Sample'].values[0]
        avg_frag_err_metrics = metrics_results['Average Fragment Error'].values[0]
        prec_error_metrics = metrics_results['Precursor Error'].values[0]
        hyperscore_metrics = metrics_results['Hyperscore'].values[0]
        consec_b_metrics = metrics_results['# Consecutive b-ions'].values[0]
        consec_y_metrics = metrics_results['# Consecutive y-ions'].values[0]
        motif_metrics = metrics_results['Motif_Score'].values[0]
        annot_per_ion_metrics = metrics_results['Average annotations/fragment '].values[0]
        frag_ions_AA_metrics = metrics_results['Average number of fragment ions per AA'].values[0]
        non_neut_per_AA_metrics = metrics_results['Number of non-neutral-loss fragment ions per AA'].values[0]
        percent_seq_cov_metrics = metrics_results['% Seq Coverage'].values[0]
        prec_err_expand = prec_error_metrics * precursor_err_mult
        consec_b_expand = consec_b_metrics * no_consec_b_ions_mult
        percent_seq_cov_expand = percent_seq_cov_metrics * seq_cov_mult
        
        multiply_list = []
        divide_list = []
 
        if consec_b_fact == 'Multiply':
            if consec_b_expand > 0:
                multiply_list.append(consec_b_expand)
            else:
                pass
        elif consec_b_fact == 'Divide':
            if consec_b_expand > 0:
                divide_list.append(consec_b_expand)
            else:
                pass
        else:
            raise ValueError('Function not specified: consec_b_fact') 

        if percent_seq_cov_fact == 'Multiply':
            if percent_seq_cov_expand > 0:
                multiply_list.append(percent_seq_cov_expand)
            else:
                pass
        elif percent_seq_cov_fact == 'Divide':
            if percent_seq_cov_expand > 0:
                divide_list.append(percent_seq_cov_expand)
            else:
                pass
        else:
            raise ValueError('Function not specified: percent_seq_cov_fact')       

        if prec_err_fact == 'Multiply':
            if prec_err_expand > 0:
                multiply_list.append(prec_err_expand)
            else:
                pass
        elif prec_err_fact == 'Divide':
            if prec_err_expand > 0:
                divide_list.append(prec_err_expand)
            else:
                pass
        else:
            raise ValueError('Function not specified: prec_err_fact')  
            
        if len(multiply_list)>1:
            numerator = np.prod(multiply_list)
        else:
            numerator = 1
        
        if len(divide_list)>1:
            denomenator = np.prod(divide_list)
        else:
            denomenator = 1
            
        final_score = numerator / denomenator
        
        metrics_results['Final score, run: 1'] = final_score
        metrics_results['Sample'] = sample
        final_score_storage.append(metrics_results)
    
    dsd_results_summary = pd.concat(final_score_storage,ignore_index=True)
    
    metrics_w_dsd_results = pd.merge(metrics_w_dsd_results,dsd_results_summary, on=['Sample','Peptide','Scan','Average Fragment Error','Precursor Error','Hyperscore',
     '# Consecutive b-ions','# Consecutive y-ions','Motif_Score','Average annotations/fragment ',
     'Average number of fragment ions per AA','Number of non-neutral-loss fragment ions per AA','% Seq Coverage'])
    
    file_path = sample_output_directory + '\\IDed_peptide_scores.csv'
    with open(file_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            metrics_w_dsd_results.to_csv(filec,index=False)
    
    return metrics_w_dsd_results
    
