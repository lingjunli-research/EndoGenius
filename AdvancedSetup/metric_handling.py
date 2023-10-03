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
    
    #dsd_path = r"C:\Users\lawashburn\Documents\DB_pep_validation\DSD_pt2_20230701\dsd_pt2_v01_20230701.csv"
    #metrics_extracted_path = r"C:\Users\lawashburn\Documents\DB_pep_validation\motif_percent_assess_2023010\metrics_extracted_motif_100_percent.csv"
    #output_directory = r"C:\Users\lawashburn\Documents\DB_pep_validation\motif_percent_assess_2023010"
    
    #dsd = pd.read_csv(dsd_path)
    metrics = metrics_extracted_path
    metrics = metrics.drop(columns=['Pyro-glu on E','Pyro-glu on Q','Oxidation','Sulfation', '# Modifications','C-termini amidation',
                                    'Max Fragment Error','Min Fragment Error','Median Fragment Error','% Fragment ions are b',
                                    '% Fragment ions are y'])
    
    metrics_w_dsd_results = metrics
    
    
    run = 1
    #frag_err_mult = dsd_run['Avg Fragment Error'].values[0]
    precursor_err_mult = 10
    no_consec_b_ions_mult = 10
    #no_consec_y_ions_mult = dsd_run['# Consecutive y-ions'].values[0]
    #annot_per_frag_mult = dsd_run['Avg annotations/fragment'].values[0]
    seq_cov_mult = 10
    #frag_ions_per_AA_mult = dsd_run['Avg frag ions/AA'].values[0]
    non_neut_frag_ions_per_AA_mult = 10
    #hyperscore_mult = dsd_run['Hyperscore'].values[0]
    #motif_score_mult = dsd_run['Motif score'].values[0]
    
    # hyperscore_fact = dsd_run['Hyperscore_factor'].values[0]
    # avg_frag_err_fact = dsd_run['Avg frag error mult'].values[0]
    consec_b_fact = 'Multiply'
    # consec_y_fact = dsd_run['Consec yions mult'].values[0]
    # #annot_per_ion_fact = dsd_run['Annotate per ions mult'].values[0]
    percent_seq_cov_fact = 'Multiply'
    # #frag_ions_AA_fact = dsd_run['Avg frag ions per AA mult'].values[0]
    non_neut_per_AA_fact = 'Divide'
    # motif_fact = dsd_run['Motif mult'].values[0]
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
        
        #avg_frag_err_expand = avg_frag_err_metrics * frag_err_mult
        prec_err_expand = prec_error_metrics * precursor_err_mult
        #hyperscore_expand = hyperscore_metrics * hyperscore_mult
        consec_b_expand = consec_b_metrics * no_consec_b_ions_mult
        #consec_y_expand = consec_y_metrics * no_consec_y_ions_mult
        # motif_expand = motif_metrics * motif_score_mult
        # annot_per_ion_expand = annot_per_ion_metrics * annot_per_frag_mult
        # frag_ions_AA_expand = frag_ions_AA_metrics * frag_ions_per_AA_mult
        #non_neut_per_AA_expand = non_neut_per_AA_metrics * non_neut_frag_ions_per_AA_mult
        percent_seq_cov_expand = percent_seq_cov_metrics * seq_cov_mult
        
        multiply_list = []
        divide_list = []
        
        # if hyperscore_fact == 'Multiply':
        #     if hyperscore_expand > 0:
        #         multiply_list.append(hyperscore_expand)
        #     else:
        #         pass
        # elif hyperscore_fact == 'Divide':
        #     if hyperscore_expand > 0:
        #         divide_list.append(hyperscore_expand)
        #     else:
        #         pass
        # else:
        #     raise ValueError('Function not specified: hyperscore_fact')
            
        # if avg_frag_err_fact == 'Multiply':
        #     if avg_frag_err_expand > 0:
        #         multiply_list.append(avg_frag_err_expand)
        #     else:
        #         pass
        # elif avg_frag_err_fact == 'Divide':
        #     if avg_frag_err_expand > 0:
        #         divide_list.append(avg_frag_err_expand)
        #     else:
        #         pass
        # else:
        #     raise ValueError('Function not specified: avg_frag_err_fact')
            
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
            
        # if consec_y_fact == 'Multiply':
        #     if consec_y_expand > 0:
        #         multiply_list.append(consec_y_expand)
        #     else:
        #         pass
        # elif consec_y_fact == 'Divide':
        #     if consec_y_expand > 0:
        #         divide_list.append(consec_y_expand)
        #     else:
        #         pass
        # else:
        #     raise ValueError('Function not specified: consec_y_fact')  
            
        # if annot_per_ion_fact == 'Multiply':
        #     if annot_per_ion_expand > 0:
        #         multiply_list.append(annot_per_ion_expand)
        #     else:
        #         pass
        # elif annot_per_ion_fact == 'Divide':
        #     if annot_per_ion_expand > 0:
        #         divide_list.append(annot_per_ion_expand)
        #     else:
        #         pass
        # else:
        #     raise ValueError('Function not specified: annot_per_ion_fact')   
            
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
            
        # if frag_ions_AA_fact == 'Multiply':
        #     if frag_ions_AA_expand > 0:
        #         multiply_list.append(frag_ions_AA_expand)
        #     else:
        #         pass
        # elif frag_ions_AA_fact == 'Divide':
        #     if frag_ions_AA_expand > 0:
        #         divide_list.append(frag_ions_AA_expand)
        #     else:
        #         pass
        # else:
        #     raise ValueError('Function not specified: frag_ions_AA_fact')     
            
        # if non_neut_per_AA_fact == 'Multiply':
        #     if non_neut_per_AA_expand > 0:
        #         multiply_list.append(non_neut_per_AA_expand)
        #     else:
        #         pass
        # elif non_neut_per_AA_fact == 'Divide':
        #     if non_neut_per_AA_expand > 0:
        #         divide_list.append(non_neut_per_AA_expand)
        #     else:
        #         pass
        # else:
        #     raise ValueError('Function not specified: non_neut_per_AA_fact')     
            
        # if motif_fact == 'Multiply':
        #     if motif_expand > 0:
        #         multiply_list.append(motif_expand)
        #     else:
        #         pass
        # elif motif_fact == 'Divide':
        #     if motif_expand > 0:
        #         divide_list.append(motif_expand)
        #     else:
        #         pass
        # else:
        #     raise ValueError('Function not specified: motif_fact')       
    
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
    
