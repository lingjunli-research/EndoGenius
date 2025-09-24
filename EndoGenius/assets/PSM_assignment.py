# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 10:41:46 2023

@author: lawashburn
"""

import pandas as pd
import csv
import numpy as np
import re
import os
import smtplib
import itertools

def PSM_assignment_execute(standard_err_percent,confident_seq_cov,max_adjacent_swapped_AA,min_motif_len,fragment_error_threshold,num_sub_AAs,db_search_parent_directory,target_path,motif_path):

    # standard_err_percent = 0.1
    # confident_seq_cov = 70
    # max_adjacent_swapped_AA = 2
    # min_motif_len = 3
    # fragment_error_threshold = 0.02
    # num_sub_AAs = 1
    
    # db_search_parent_directory = r"C:\Users\lawashburn\Documents\DB_pep_validation\DB_search_w_mod_20230629"
    
    # target_path = r"C:\Users\lawashburn\Documents\DBpep_v2\finale\Reference_DB\target_list.csv"
    target = pd.read_csv(target_path)
    target_list =  target['Sequence'].values.tolist()
    
    # motif_path = r"C:\Users\lawashburn\Documents\DBpep_v2\Validation_w_Kellen_Motif_Nhu_Raw_Files\PEAKS_oursoftware_compare_brain_only\DB_search_optimize\03_precursor_AMM_optimize_v2\2021_0817_Brain_1\motif_db_kellen_tina__lauren_20230621_v03.csv"
    motif_db = pd.read_csv(motif_path)
    
    motif_list = motif_db['Sequence'].values.tolist()
    
    def compare(string1, string2, no_match_c=' ', match_c='|'):
        if len(string2) < len(string1):
            string1, string2 = string2, string1
        result = ''
        n_diff = 0
        for c1, c2 in zip(string1, string2):
            if c1 == c2:
                result += match_c
            else:
                result += no_match_c
                n_diff += 1
        delta = len(string2) - len(string1)
        result += delta * no_match_c
        n_diff += delta
        return (result, n_diff)
    def dif_compare(string1, string2):
        result, n_diff = compare(string1, string2, no_match_c='_')
        return(n_diff)
    def get_dir_names_with_strings_list(str_list): #definition for finding a file containing a string in filename in specified directory
        full_list = os.listdir(db_search_parent_directory)
        final_list = [nm for ps in str_list for nm in full_list if ps in nm]
        return final_list
    
    query = '' #search for ion list pertaining to the sequence
    parent_dir_list = (get_dir_names_with_strings_list([query])) #search for the file based on query
    
    for directory in parent_dir_list:
        all_correlation_results_path = db_search_parent_directory + '\\' + directory +'\\all_correlation_results.csv'
        all_correlation_results = pd.read_csv(all_correlation_results_path)
        
        all_correlation_results['Sequence with mod'] = all_correlation_results['Sequence']
        all_correlation_results['Sequence'] = all_correlation_results['Sequence'].str.replace(r"\([^()]*\)", "", regex=True)
        
        peptide_report_output = db_search_parent_directory + '\\' + directory
    
        ##Begin PSM Assignment##
        max_seq_cov = all_correlation_results['Sequence coverage'].max()
        max_thresh = max_seq_cov * (1-standard_err_percent)
        
        standard_err_subset = all_correlation_results[all_correlation_results['Sequence coverage'] >= max_thresh]
        standard_err = standard_err_subset['Correlation value'].std()
        
        all_correlation_results['count'] = all_correlation_results.groupby('Scan')['Scan'].transform(len)
        
        final_psm = all_correlation_results[all_correlation_results['count'] == 1]
        final_psm['Step assigned'] = 'Only peptide match for scan'
        
        psm_candidate = all_correlation_results[all_correlation_results['count'] > 1]
        
        sequence_list = psm_candidate['Sequence'].values.tolist()
        
        seq_log = []
        motif_log = []
        status_log = []
        
        for a in sequence_list:
            if a not in seq_log:
                for b in motif_list:
                    if b in a:
                        seq_log.append(a)
                        motif_log.append(b)
                        status_log.append(True)
                    else:
                        pass
            else:
                pass
            
        for c in sequence_list:
            if c not in seq_log:
                seq_log.append(c)
                motif_log.append(np.nan)
                status_log.append(False)
                    
        motif_eval = pd.DataFrame()
        motif_eval['Sequence'] = seq_log
        motif_eval['Motif'] = motif_log
        motif_eval['Motif status'] = status_log
        
        psm_candidate = psm_candidate.merge(motif_eval,on='Sequence')
        
        scan_candidate = psm_candidate.drop_duplicates(subset='Scan')
        scan_candidate_list = scan_candidate['Scan'].values.tolist()
        
        future_storage = []
        
        for scan in scan_candidate_list:
            final_psm = final_psm.drop_duplicates()
            scan_filtered_psm_candidate = psm_candidate[psm_candidate['Scan'] == scan]
            filter_A = scan_filtered_psm_candidate[scan_filtered_psm_candidate['Motif status'] == True]
            if len(filter_A) == 1:
                filter_A['Step assigned'] = 'filter A(+)'
                final_psm = pd.concat([final_psm,filter_A])
            elif len(filter_A) > 1:
                criteria_AB = filter_A
                criteria_AB = criteria_AB.copy()
                criteria_AB['Seq Length']  = criteria_AB['Sequence'].str.len()
                criteria_AB = criteria_AB.copy()
                criteria_AB['Motif Length']  = criteria_AB['Motif'].str.len()
                criteria_AB = criteria_AB.copy()
                criteria_AB['Motif:Seq Ratio'] = criteria_AB['Motif Length'] * np.sqrt(criteria_AB['Seq Length'])
                value_AB = criteria_AB['Motif:Seq Ratio'].max()
                criteria_AB = criteria_AB[criteria_AB['Motif:Seq Ratio'] == value_AB]
                if len(criteria_AB) == 1:
                    sequence_AB = criteria_AB['Sequence'].values[0]
                    filter_AB = filter_A[filter_A['Sequence'] == sequence_AB]
                    filter_AB['Step assigned'] = 'filter A(++) B(+)'
                    final_psm = pd.concat([final_psm,filter_AB])
                if len(criteria_AB) > 1:
                    entry_storage = []
                    
                    for k in range(0,len(criteria_AB)):
                        entry = criteria_AB.iloc[[k]]
                        sequence_to_check = entry['Sequence'].values[0]
                        
                        motif_count_storage = []
                        for l in motif_list:
                            if l in sequence_to_check:
                                motif_count_storage.append(l)
                        entry = entry.copy()
                        entry['# Matching motifs'] = len(motif_count_storage)
                        entry_storage.append(entry)
                    
                    motif_w_counts = pd.concat(entry_storage,ignore_index=True) 
                    
                    value_ABC = motif_w_counts['# Matching motifs'].max()
                    filter_ABC = motif_w_counts[motif_w_counts['# Matching motifs'] == value_ABC]
                    if len(filter_ABC) == 1: 
                        sequence_ABC = filter_ABC['Sequence'].values[0]
                        psm_entry = criteria_AB[criteria_AB['Sequence'] == sequence_ABC] #this one was previously criteria A not criteria AB
                        psm_entry['Step assigned'] = 'filter A(++) B(++) C(+)'
                        final_psm = pd.concat([final_psm,psm_entry])
                    
                    sequence_storage2 = []
                    sequence_motif_coverage = []
                    if len(filter_ABC) > 1:
                        for k in range(0,len(filter_ABC)):
                            entry = filter_ABC.iloc[[k]]
                            sequence_to_check = entry['Sequence'].values[0]
                            
                            motif_index_storage = []
                            for l in motif_list:
                                if l in sequence_to_check:
                                    for match in re.finditer(l, sequence_to_check):
                                        start = match.start()
                                        end = match.end()
                                        for index in range(start,end):
                                            if index not in motif_index_storage:
                                                motif_index_storage.append(index)
                            
                            seq_cov_motif = (len(motif_index_storage))/(len(sequence_to_check))
                            sequence_motif_coverage.append(seq_cov_motif)
                            sequence_storage2.append(sequence_to_check)
                    
                    sequence_coverage_motif_df = pd.DataFrame()
                    sequence_coverage_motif_df['Sequence'] = sequence_storage2
                    sequence_coverage_motif_df['Sequence/Motif Coverage'] = sequence_motif_coverage
                    sequence_coverage_motif_df = sequence_coverage_motif_df.drop_duplicates()
                    
                    value_ABCD = sequence_coverage_motif_df['Sequence/Motif Coverage'].max()
                    criteria_ABCD = sequence_coverage_motif_df[sequence_coverage_motif_df['Sequence/Motif Coverage'] == value_ABCD]
        
                    if len(criteria_ABCD) == 1:
                        sequence_ABCD = criteria_ABCD['Sequence'].values[0]
                        psm_entry = filter_ABC[filter_ABC['Sequence'] == sequence_ABCD]  #this one was previously filter A not filter ABC
                        psm_entry['Step assigned'] = 'filter A(++) B(++) C(++) D(+)'
                        final_psm = pd.concat([final_psm,psm_entry])
                    
                    elif len(criteria_ABCD) > 1:     
                        filter_list = criteria_ABCD['Sequence'].values.tolist()
                        criteria_ABCDE = filter_ABC[filter_ABC['Sequence'].isin(filter_list)]
                        filter_ABCDE = criteria_ABCDE[criteria_ABCDE['Sequence coverage'] >= confident_seq_cov]
                        if len(filter_ABCDE) == 1:
                            filter_ABCDE['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(+)'
                            final_psm = pd.concat([final_psm,filter_ABCDE])
                        
                        elif len(filter_ABCDE) > 1:
                            value_ABCDEG = filter_ABCDE['Correlation value'].max()
                            filter_ABCDEG = filter_ABCDE[filter_ABCDE['Correlation value'] == value_ABCDEG]
                            filter_ABCDEG = filter_ABCDEG.drop_duplicates(subset='Sequence')
                            if len(filter_ABCDEG) == 1:
                                filter_ABCDEG['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(++) G(+)'
                                final_psm = pd.concat([final_psm,filter_ABCDEG])
                            if len(filter_ABCDEG) > 1:
                                value_ABCDEGI = filter_ABCDEG['Sequence coverage'].max()
                                filter_ABCDEGI = filter_ABCDEG[filter_ABCDEG['Sequence coverage'] == value_ABCDEGI]
                                if len(filter_ABCDEGI) == 1:
                                    filter_ABCDEGI['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(++) G(++) I(+)'
                                    final_psm = pd.concat([final_psm,filter_ABCDEGI])
                                if len(filter_ABCDEGI) > 1:
    
                                    filter_ABCDEGI['Seq_Len']  = filter_ABCDEGI['Sequence'].str.len()
                                    max_seq_len = filter_ABCDEGI['Seq_Len'].max()
                                    criteria_ABCDEGIH = filter_ABCDEGI[filter_ABCDEGI['Seq_Len'] == max_seq_len]
                                    if len(filter_ABCDEGI) == len(criteria_ABCDEGIH): #filter for peptides of same length
                                        peptide_list_same_len = criteria_ABCDEGIH['Sequence'].values.tolist()
                                        if len(peptide_list_same_len) == 2:
                                            seq1 = peptide_list_same_len[0]
                                            seq2 = peptide_list_same_len[1]
                                            differing_indexes = indexes = [i for i in range(len(seq1)) if seq1[i] != seq2[i]]
                                            s = pd.Series(differing_indexes)
                                            sgc = s.groupby(s.diff().ne(1).cumsum()).transform('count')
                                            result = s[sgc == sgc.max()].tolist()
                                            if len(result) <= max_adjacent_swapped_AA:
                                                filter_ABCDEGIH = criteria_ABCDEGIH[criteria_ABCDEGIH['Sequence'].isin(peptide_list_same_len)]
                                                filter_ABCDEGIH['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(++) G(++) I(++) H(+)'
                                                final_psm = pd.concat([final_psm,filter_ABCDEGIH])
                                            if len(result) > max_adjacent_swapped_AA:
                                                #raise ValueError('Check this') #no issue, validated 7/8/23
                                                psm_entry = filter_ABCDEGI.sample(n=1, random_state=1)
                                                psm_entry['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(++) G(++) I(++) H(-?)'
                                                final_psm = pd.concat([final_psm,psm_entry])
                                        if len(peptide_list_same_len) != 2:
                                            #raise ValueError('Check this') #no impact, validated 7/10/23-2
                                            psm_entry = filter_ABCDEGI.sample(n=1, random_state=1)
                                            psm_entry['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(++) G(++) I(++) H(-?)'
                                            final_psm = pd.concat([final_psm,psm_entry])
                                    if len(filter_ABCDEGI) != len(criteria_ABCDEGIH): 
                                         #raise ValueError('Check this') #no impact, validated 7/10/23-2
                                         psm_entry = filter_ABCDEGI.sample(n=1, random_state=1)
                                         psm_entry['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(++) G(++) I(++) H(-?)'
                                         final_psm = pd.concat([final_psm,psm_entry])
                                   
                        
                        elif len(filter_ABCDE) == 0:
                            filter_list = criteria_ABCD['Sequence'].values.tolist()
                            criteria_ABCDEGF = filter_ABC[filter_ABC['Sequence'].isin(filter_list)]
                            criteria_ABCDEGF = criteria_ABCDEGF.drop_duplicates(subset='Sequence')
                            filter_ABCDEGF = criteria_ABCDEGF[criteria_ABCDEGF['Correlation value'] > 0]
                            
                            if len(filter_ABCDEGF) == 1:
                                filter_ABCDEGF['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(-) G(-) F(+)'
                                final_psm = pd.concat([final_psm,filter_ABCDEGF])
                            if len(filter_ABCDEGF) > 1: 
                                value_ABCDEGFG = filter_ABCDEGF['Correlation value'].max()
                                filter_ABCDEGFG = filter_ABCDEGF[filter_ABCDEGF['Correlation value'] == value_ABCDEGFG]
                                if len(filter_ABCDEGFG) == 1:
                                    filter_ABCDEGFG['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(-) G(++) F(++) G(+)'
                                    final_psm = pd.concat([final_psm,filter_ABCDEGFG])
                                if len(filter_ABCDEGFG) > 1:
                                    filter_ABCDEGFG['Seq_Len']  = filter_ABCDEGFG['Sequence'].str.len()
                                    max_seq_len = filter_ABCDEGFG['Seq_Len'].max()
                                    criteria_ABCDEGFGH = filter_ABCDEGFG[filter_ABCDEGFG['Seq_Len'] == max_seq_len]
                                    if len(filter_ABCDEGFG) == len(criteria_ABCDEGFGH): #filter for peptides of same length
                                        peptide_list_same_len = criteria_ABCDEGFGH['Sequence'].values.tolist()
                                        if len(peptide_list_same_len) == 2:
                                            seq1 = peptide_list_same_len[0]
                                            seq2 = peptide_list_same_len[1]
                                            differing_indexes = indexes = [i for i in range(len(seq1)) if seq1[i] != seq2[i]]
                                            s = pd.Series(differing_indexes)
                                            sgc = s.groupby(s.diff().ne(1).cumsum()).transform('count')
                                            result = s[sgc == sgc.max()].tolist()
                                            if len(result) <= max_adjacent_swapped_AA:
                                                filter_ABCDEGFGH = criteria_ABCDEGFGH[criteria_ABCDEGFGH['Sequence'].isin(peptide_list_same_len)]
                                                filter_ABCDEGFGH['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(-) G(++) F(++) G(++) H(+)'
                                                final_psm = pd.concat([final_psm,filter_ABCDEGFGH])
                                            if len(result) > max_adjacent_swapped_AA:
                                                
                                                scan_store = []
                                                seq_store = []
                                                len_store = []
                                                
                                                for line in range(0,len(filter_ABCDEGFG)):
                                                    criteria_ABCDEGFGJ = filter_ABCDEGFG.iloc[[line]]
                                                    scan_oi = criteria_ABCDEGFGJ['Scan'].values[0]
                                                    seq_oi = criteria_ABCDEGFGJ['Sequence'].values[0]
                                                    seq_mod = criteria_ABCDEGFGJ['Sequence with mod'].values[0]
                                                    seq_mod = seq_mod.replace('(Glu->pyro-Glu)','(pyroGlu)')
                                                    seq_mod = seq_mod.replace('(Gln->pyro-Glu)','(pyroGlu)')
                                                    
                                                    fragment_report_path = peptide_report_output + '\\fragment_matches\\' + seq_mod + '_' + str(scan_oi) + '_fragment_report.csv'
                                                    fragment_report = pd.read_csv(fragment_report_path)
                                                    fragment_report = fragment_report[fragment_report['Fragment error (Da)'] >= fragment_error_threshold]
                                                    
                                                    
                                                    fragment_report_len = len(fragment_report)
                                                    
                                                    scan_store.append(scan_oi)
                                                    seq_store.append(seq_oi)
                                                    len_store.append(fragment_report_len)
                                                
                                                scan_seq_len = pd.DataFrame()
                                                scan_seq_len['Sequence'] = seq_store
                                                scan_seq_len['Scan'] = scan_store
                                                scan_seq_len['Length'] = len_store
                                                
                                                value_ABCDEGFGJ = scan_seq_len['Length'].max()
                                                criteria_ABCDEGFGJ2 = scan_seq_len[scan_seq_len['Length'] == value_ABCDEGFGJ]
                                                criteria_str_ABCDEGFGJ2 = criteria_ABCDEGFGJ2['Sequence'].values.tolist()
                                                
                                                filter_ABCDEGFGJ = filter_ABCDEGFG[filter_ABCDEGFG['Sequence'].isin(criteria_str_ABCDEGFGJ2)]
                                                if len(filter_ABCDEGFGJ) == 1:
                                                    filter_ABCDEGFGJ['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(-) G(++) F(++) G(++) H(-) J(+)'
                                                    final_psm = pd.concat([final_psm,filter_ABCDEGFGJ])
                                                if len(filter_ABCDEGFGJ) > 1:
                                                    #raise ValueError('Check this') #no impact, validated 7/10/23-2
                                                    psm_entry = filter_ABCDEGFGJ.sample(n=1, random_state=1)
                                                    psm_entry['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(-) G(++) F(++) G(++) H(-) J(++?)'
                                                    final_psm = pd.concat([final_psm,psm_entry])
        
                                        else:
                                            #raise ValueError('Check this') #no impact, validated 7/10/23-2
                                            psm_entry = filter_ABCDEGFG.sample(n=1, random_state=1)
                                            psm_entry['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(-) G(++) F(++) G(++) H(-?)'
                                            final_psm = pd.concat([final_psm,psm_entry])
                                        
                                    if len(filter_ABCDEGFG) != len(criteria_ABCDEGFGH): #filter for peptides of same length  
                                        #raise ValueError('Check this') #no impact, validated 7/10/23-2
                                        psm_entry = filter_ABCDEGFG.sample(n=1, random_state=1)
                                        psm_entry['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(-) G(++) F(++) G(++) H(-?)'
                                        final_psm = pd.concat([final_psm,psm_entry])
                            if len(filter_ABCDEGF) == 0:   
                                value_ABCDEGFI = filter_ABC['Sequence coverage'].max()
                                filter_ABCDEGFI = filter_ABC[filter_ABC['Sequence coverage'] == value_ABCDEGFI]
                                if len(filter_ABCDEGFI) == 1:
                                    filter_ABCDEGFI['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(-) G(++) F(++) I(+)'
                                    final_psm = pd.concat([final_psm,filter_ABCDEGFI])
                                if len(filter_ABCDEGFI) > 1:
                                    scan_store = []
                                    seq_store = []
                                    len_store = []
                                    
                                    for line in range(0,len(filter_ABCDEGFI)):
                                        criteria_ABCDEGFIJ = filter_ABCDEGFI.iloc[[line]]
                                        scan_oi = criteria_ABCDEGFIJ['Scan'].values[0]
                                        seq_oi = criteria_ABCDEGFIJ['Sequence'].values[0]
                                        
                                        seq_mod = criteria_ABCDEGFIJ['Sequence with mod'].values[0]
                                        seq_mod = seq_mod.replace('(Glu->pyro-Glu)','(pyroGlu)')
                                        seq_mod = seq_mod.replace('(Gln->pyro-Glu)','(pyroGlu)')
                                        
                                        fragment_report_path = peptide_report_output + '\\fragment_matches\\' + seq_mod + '_' + str(scan_oi) + '_fragment_report.csv'
                                        fragment_report = pd.read_csv(fragment_report_path)
                                        fragment_report = fragment_report[fragment_report['Fragment error (Da)'] >= fragment_error_threshold]
                                        
                                        
                                        fragment_report_len = len(fragment_report)
                                        
                                        scan_store.append(scan_oi)
                                        seq_store.append(seq_oi)
                                        len_store.append(fragment_report_len)
                                    
                                    scan_seq_len = pd.DataFrame()
                                    scan_seq_len['Sequence'] = seq_store
                                    scan_seq_len['Scan'] = scan_store
                                    scan_seq_len['Length'] = len_store
                                    
                                    value_ABCDEGFIJ = scan_seq_len['Length'].max()
                                    criteria_ABCDEGFIJ2 = scan_seq_len[scan_seq_len['Length'] == value_ABCDEGFIJ]
                                    criteria_str_ABCDEGFIJ2 = criteria_ABCDEGFIJ2['Sequence'].values.tolist()
                                    
                                    filter_ABCDEGFIJ = filter_ABCDEGFI[filter_ABCDEGFI['Sequence'].isin(criteria_str_ABCDEGFIJ2)]
                                    filter_ABCDEGFIJ = filter_ABCDEGFIJ.drop_duplicates()
                                    if len(filter_ABCDEGFIJ) == 1:
                                        filter_ABCDEGFIJ['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(-) G(++) F(++) I(++) J(++)'
                                        final_psm = pd.concat([final_psm,filter_ABCDEGFIJ])
                                    if len(filter_ABCDEGFIJ) > 1:
                                        filter_ABCDEGFIJ['Seq_Len']  = filter_ABCDEGFIJ['Sequence'].str.len()
                                        max_seq_len = filter_ABCDEGFIJ['Seq_Len'].max()
                                        criteria_ABCDEGFIJH = filter_ABCDEGFIJ[filter_ABCDEGFIJ['Seq_Len'] == max_seq_len]
                                        if len(filter_ABCDEGFIJ) == len(criteria_ABCDEGFIJH): #filter for peptides of same length
                                            peptide_list_same_len = criteria_ABCDEGFIJH['Sequence'].values.tolist()
                                            if len(peptide_list_same_len) == 2:
                                                seq1 = peptide_list_same_len[0]
                                                seq2 = peptide_list_same_len[1]
                                                differing_indexes = indexes = [i for i in range(len(seq1)) if seq1[i] != seq2[i]]
                                                s = pd.Series(differing_indexes)
                                                sgc = s.groupby(s.diff().ne(1).cumsum()).transform('count')
                                                result = s[sgc == sgc.max()].tolist()
                                                if len(result) <= max_adjacent_swapped_AA:
                                                    filter_ABCDEGFIJH = criteria_ABCDEGFIJH[criteria_ABCDEGFIJH['Sequence'].isin(peptide_list_same_len)]
                                                    filter_ABCDEGFIJH['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(-) G(++) F(++) I(++) J(++) H(+)'
                                                    final_psm = pd.concat([final_psm,filter_ABCDEGFIJH])
                                                if len(result) > max_adjacent_swapped_AA:
                                                    motif_storage = []
                                                    seq_storage = []
                                                    position_storage = []
    
                                                    for ind in range(0,len(filter_ABCDEGFIJ)):
                                                        ind_filter = filter_ABCDEGFIJ.iloc[[ind]]
                                                        motif = ind_filter['Motif'].values[0]
                                                        seq = ind_filter['Sequence'].values[0]
                                                        
                                                        if seq.startswith((motif)):
                                                            position_storage.append(True)
                                                            motif_storage.append(motif)
                                                            seq_storage.append(seq)
                                                        elif seq.endswith((motif)):
                                                            position_storage.append(True)
                                                            motif_storage.append(motif)
                                                            seq_storage.append(seq)   
                                                        else:
                                                            position_storage.append(False)
                                                            motif_storage.append(motif)
                                                            seq_storage.append(seq)
                                                        
                                                    position_df = pd.DataFrame()
                                                    position_df['Sequence'] = seq_storage
                                                    position_df['Motif'] = motif_storage
                                                    position_df['Position'] = position_storage
                                                    
                                                    filter_ABCDEGFIJM = position_df[position_df['Position'] == True]
                                                    
                                                    if len(filter_ABCDEGFIJM) == 0:
                                                        #raise ValueError('Check this') #no impact, validated 7/10/23-2
                                                        psm_entry = filter_ABCDEGFIJ.sample(n=1, random_state=1)
                                                        psm_entry['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(-) G(++) F(++) I(++) J(++) H(-) M(-)'
                                                        final_psm = pd.concat([final_psm,psm_entry])
                                                        
                                                    if len(filter_ABCDEGFIJM) == 1:
                                                        value_ABCDEGFIJM = filter_ABCDEGFIJM['Sequence'].values[0]
                                                        filter_ABCDEGFIJM2 = filter_ABCDEGFIJ[filter_ABCDEGFIJ['Sequence'] == value_ABCDEGFIJM]
                                                        filter_ABCDEGFIJM2['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(-) G(++) F(++) I(++) J(++) H(-) M(+)'
                                                        final_psm = pd.concat([final_psm,filter_ABCDEGFIJM2])
                                                    
                                                    if len(filter_ABCDEGFIJM) > 1:  
                                                        #raise ValueError('Check this') #no impact, validated 7/10/23-2
                                                        value_ABCDEGFIJM = filter_ABCDEGFIJM['Sequence'].values.tolist()
                                                        filter_ABCDEGFIJM2 = filter_ABCDEGFIJM[filter_ABCDEGFIJM['Sequence'].isin(value_ABCDEGFIJM)]
                                                        psm_entry = filter_ABCDEGFIJM2.sample(n=1, random_state=1)
                                                        psm_entry['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(-) G(++) F(++) I(++) J(++) H(-) M(++?)'
                                                        final_psm = pd.concat([final_psm,psm_entry])
                                                   
                                            if len(peptide_list_same_len) != 2:
                                                filter_ABCDEGFIJN = filter_ABCDEGFIJ.drop_duplicates(subset='Sequence')
                                                if len(filter_ABCDEGFIJN) == 1:
                                                    filter_ABCDEGFIJN['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(-) G(++) F(++) I(++) J(+) N(+)'
                                                    final_psm = pd.concat([final_psm,filter_ABCDEGFIJN])
                                                if len(filter_ABCDEGFIJN) > 1:
                                                    #raise ValueError('Check this') #no impact, validated 7/8/23
                                                    psm_entry = filter_ABCDEGFIJN.sample(n=1, random_state=1)
                                                    psm_entry['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(-) G(++) F(++) I(++) J(++) N(++?)'
                                                    final_psm = pd.concat([final_psm,psm_entry])
                                        if len(filter_ABCDEGFIJ) != len(criteria_ABCDEGFIJH):
                                            #raise ValueError('Check this') #no impact, validated 7/10/23
                                            psm_entry = filter_ABCDEGFIJ.sample(n=1, random_state=1)
                                            psm_entry['Step assigned'] = 'filter A(++) B(++) C(++) D(++) E(-) G(++) F(++) I(++) J(++) H(-)'
                                            final_psm = pd.concat([final_psm,psm_entry])
        
            else:
                # 
                filter_AE = scan_filtered_psm_candidate[scan_filtered_psm_candidate['Sequence coverage'] >= confident_seq_cov]
                
                if len(filter_AE) == 1:
                    filter_AE['Step assigned'] = 'filter A(-) E(+)'
                    final_psm = pd.concat([final_psm,filter_AE])
                    
                if len(filter_AE)>1:
                    value_AEI = filter_AE['Sequence coverage'].max()
                    filter_AEI = filter_AE[filter_AE['Sequence coverage'] == value_AEI]
                    if len(filter_AEI) == 1:
                        filter_AEI['Step assigned'] = 'filter A(-) E(+) I(+)'
                        final_psm = pd.concat([final_psm,filter_AEI])
                    if len(filter_AEI) > 1: 
                        filter_AEIK = filter_AEI
                        filter_AEIK['Correlation value (rounded)'] = filter_AEIK['Correlation value'].round(0)
                        value_AEIK = filter_AEIK['Correlation value (rounded)'].max()
                        filter_AEIK = filter_AEIK[filter_AEIK['Correlation value (rounded)'] == value_AEIK]
                        filter_AEIK = filter_AEIK.drop(columns = 'Correlation value (rounded)')
                        if len(filter_AEIK) == 1:
                            
                            filter_AEIK['Step assigned'] = 'filter A(-) E(++) I(++) K(+)'
                            final_psm = pd.concat([final_psm,filter_AEIK])
                        if len(filter_AEIK) > 1:
                            filter_AEIK['Seq_Len']  = filter_AEIK['Sequence'].str.len()
                            max_seq_len = filter_AEIK['Seq_Len'].max()
                            criteria_AEIKH = filter_AEIK[filter_AEIK['Seq_Len'] == max_seq_len]
                            if len(filter_AEIK) == len(criteria_AEIKH): #filter for peptides of same length
                                peptide_list_same_len = criteria_AEIKH['Sequence'].values.tolist()
                                if len(peptide_list_same_len) == 2:
                                    seq1 = peptide_list_same_len[0]
                                    seq2 = peptide_list_same_len[1]
                                    differing_indexes = indexes = [i for i in range(len(seq1)) if seq1[i] != seq2[i]]
                                    s = pd.Series(differing_indexes)
                                    sgc = s.groupby(s.diff().ne(1).cumsum()).transform('count')
                                    result = s[sgc == sgc.max()].tolist()
                                    if len(result) <= max_adjacent_swapped_AA:
                                        filter_AEIKH = criteria_AEIKH[criteria_AEIKH['Sequence'].isin(peptide_list_same_len)]
                                        filter_AEIKH['Step assigned'] = 'filter A(-)E(++)I(++)K(++)H(+)'
                                        final_psm = pd.concat([final_psm,filter_AEIKH])
                                    if len(result) > max_adjacent_swapped_AA:
                                        #raise ValueError('Check this') #no impact, validated 7/10/23-2
                                        psm_entry = filter_AEIK.sample(n=1, random_state=1)
                                        psm_entry['Step assigned'] = 'filter A(-)E(++)I(++)K(++)H(-?)'
                                        final_psm = pd.concat([final_psm,psm_entry])
                                else:
                                    #raise ValueError('Check this') #no impact, validated 7/10/23-2
                                    psm_entry = filter_AEIK.sample(n=1, random_state=1)
                                    psm_entry['Step assigned'] = 'filter A(-)E(++)I(++)K(++)H(-?)'
                                    final_psm = pd.concat([final_psm,psm_entry])
                                
                            if len(filter_AEIK) != len(criteria_AEIKH): #filter for peptides of same length  
                                #raise ValueError('Check this') #no impact, validated 7/10/23
                                psm_entry = filter_AEIK.sample(n=1, random_state=1)
                                psm_entry['Step assigned'] = 'filter A(-)E(++)I(++)K(++)H(-?)'
                                final_psm = pd.concat([final_psm,psm_entry])
        
                if len(filter_AE)==0:
                    value_AEI = scan_filtered_psm_candidate['Sequence coverage'].max()
                    filter_AEI = scan_filtered_psm_candidate[scan_filtered_psm_candidate['Sequence coverage'] == value_AEI]
                    if len(filter_AEI) == 1:
                        filter_AEI['Step assigned'] = 'filter A(-) E(-) I(+)'
                        final_psm = pd.concat([final_psm,filter_AEI])
                    if len(filter_AEI) > 1:
                        filter_AEIF = filter_AEI[filter_AEI['Correlation value'] > 0]
                        if len(filter_AEIF) == 1:
                            filter_AEIF['Step assigned'] = 'filter A(-) E(-) I(++) F(+)'
                            final_psm = pd.concat([final_psm,filter_AEIF])
                        if len(filter_AEIF) > 1:
                            scan_store = []
                            seq_store = []
                            len_store = []
                            
                            for line in range(0,len(filter_AEIF)):
                                criteria_AEIFJ = filter_AEIF.iloc[[line]]
                                scan_oi = criteria_AEIFJ['Scan'].values[0]
                                seq_oi = criteria_AEIFJ['Sequence'].values[0]
                                
                                seq_mod = criteria_AEIFJ['Sequence with mod'].values[0]
                                seq_mod = seq_mod.replace('(Glu->pyro-Glu)','(pyroGlu)')
                                seq_mod = seq_mod.replace('(Gln->pyro-Glu)','(pyroGlu)')
                                
                                fragment_report_path = peptide_report_output + '\\fragment_matches\\' + seq_mod + '_' + str(scan_oi) + '_fragment_report.csv'
                                fragment_report = pd.read_csv(fragment_report_path)
                                fragment_report = fragment_report[fragment_report['Fragment error (Da)'] >= fragment_error_threshold]
                                
                                
                                fragment_report_len = len(fragment_report)
                                
                                scan_store.append(scan_oi)
                                seq_store.append(seq_oi)
                                len_store.append(fragment_report_len)
                            
                            scan_seq_len = pd.DataFrame()
                            scan_seq_len['Sequence'] = seq_store
                            scan_seq_len['Scan'] = scan_store
                            scan_seq_len['Length'] = len_store
                            
                            value_AEIFJ = scan_seq_len['Length'].max()
                            criteria_AEIFJ2 = scan_seq_len[scan_seq_len['Length'] == value_AEIFJ]
                            criteria_str_AEIFJ2 = criteria_AEIFJ2['Sequence'].values.tolist()
                            
                            filter_AEIFJ = filter_AEIF[filter_AEIF['Sequence'].isin(criteria_str_AEIFJ2)]
                            if len(filter_AEIFJ) == 1:
                                filter_AEIFJ['Step assigned'] = 'filter A(-) E(-) I(++) F(++) J(+)'
                                final_psm = pd.concat([final_psm,filter_AEIFJ])
                            if len(filter_AEIFJ) > 1:
                                value_AEIFJG = filter_AEIFJ['Correlation value'].max()
                                filter_AEIFJG = filter_AEIFJ[filter_AEIFJ['Correlation value'] == value_AEIFJG]
                                if len(filter_AEIFJG) == 1:
                                    filter_AEIFJG['Step assigned'] = 'filter A(-) E(-) I(++) F(++) J(++) G(+)'
                                    final_psm = pd.concat([final_psm,filter_AEIFJG])
                                if len(filter_AEIFJG) > 1:     
                                    #raise ValueError('Check this') #needs to be updated
                                    psm_entry = filter_AEIFJG.sample(n=1, random_state=1)
                                    psm_entry['Step assigned'] = 'filter A(-) E(-) I(++) F(++) J(++) G(++?)'
                                    final_psm = pd.concat([final_psm,psm_entry])
                        if len(filter_AEIF) == 0:
                            pep_storage = []
                            err_storage = []
                            for a in range(0,len(filter_AEI)):
                                row_OI = filter_AEI.iloc[[a]]
                                seq = row_OI['Sequence'].values[0]
                                seq_mod = row_OI['Sequence with mod'].values[0]
                                seq_mod = seq_mod.replace('(Glu->pyro-Glu)','(pyroGlu)')
                                seq_mod = seq_mod.replace('(Gln->pyro-Glu)','(pyroGlu)')
                                scan = row_OI['Scan'].values[0]
                                
                                peptide_report_path = peptide_report_output + '\\fragment_matches\\' + seq_mod + '_' + str(scan) + '_fragment_report.csv'         
                                peptide_report = pd.read_csv(peptide_report_path)
                                
                                peptide_report = peptide_report[peptide_report['Fragment error (Da)'] <= fragment_error_threshold]
                                if len(peptide_report)>0:
                                    avg_fragment_error = peptide_report['Fragment error (Da)'].mean()
                                    err_storage.append(avg_fragment_error)
                                    pep_storage.append(seq)
                                else:
                                    err_storage.append(1000000)
                                    pep_storage.append(seq)
                            
                            avg_err_df = pd.DataFrame()
                            avg_err_df['Sequence'] = pep_storage
                            avg_err_df['Avg Error'] = err_storage
                            
                            min_err = avg_err_df['Avg Error'].min()
                            if min_err<=fragment_error_threshold:
                                avg_err_df = avg_err_df[avg_err_df['Avg Error'] <= min_err]
                                if len(avg_err_df)==1:
                                    winning_seq = avg_err_df['Sequence'].values[0]
                                    filter_AEIL = filter_AEI[filter_AEI['Sequence'] == winning_seq]
                                    filter_AEIL['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(+)'
                                    final_psm = pd.concat([final_psm,filter_AEIL])
                                if len(avg_err_df)>1:
                                    winning_seqs = avg_err_df['Sequence'].values.tolist()
                                    filter_AEIL = filter_AEI[filter_AEI['Sequence'].isin(winning_seqs)]
                                    filter_AEILN = filter_AEIL.drop_duplicates(subset='Sequence')
                                    if len(filter_AEILN) == 1:
                                        filter_AEILN['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(++) (N+)'
                                        final_psm = pd.concat([final_psm,filter_AEILN])
                                    if len(filter_AEILN) > 1:
                                        filter_AEILN['Seq_Len']  = filter_AEILN['Sequence'].str.len()
                                        max_seq_len = filter_AEILN['Seq_Len'].max()
                                        criteria_AEILNH = filter_AEILN[filter_AEILN['Seq_Len'] == max_seq_len]
                                        if len(filter_AEILN) == len(criteria_AEILNH): #filter for peptides of same length
                                            peptide_list_same_len = criteria_AEILNH['Sequence'].values.tolist()
                                            if len(peptide_list_same_len) == 2:
                                                seq1 = peptide_list_same_len[0]
                                                seq2 = peptide_list_same_len[1]
                                                differing_indexes = indexes = [i for i in range(len(seq1)) if seq1[i] != seq2[i]]
                                                s = pd.Series(differing_indexes)
                                                sgc = s.groupby(s.diff().ne(1).cumsum()).transform('count')
                                                result = s[sgc == sgc.max()].tolist()
                                                if len(result) <= max_adjacent_swapped_AA:
                                                    filter_AEILNH = criteria_AEILNH[criteria_AEILNH['Sequence'].isin(peptide_list_same_len)]
                                                    filter_AEILNH['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(++) (N++) H(+)'
                                                    final_psm = pd.concat([final_psm,filter_AEILNH])
                                                if len(result) > max_adjacent_swapped_AA:
                                                        #raise ValueError('Check this')
                                                        psm_entry = filter_AEILN.sample(n=1, random_state=1)
                                                        psm_entry['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(++) (N++) H(-)'
                                                        final_psm = pd.concat([final_psm,psm_entry])
                                            if len(peptide_list_same_len) != 2:
                                                scan_store = []
                                                seq_store = []
                                                min_val_store = []
                                                for line in range(0,len(filter_AEI)):
                                                    criteria_AEIO = filter_AEI.iloc[[line]]
                                                    scan_oi = criteria_AEIO['Scan'].values[0]
                                                    seq_oi = criteria_AEIO['Sequence'].values[0]
                                                    
                                                    seq_mod = criteria_AEIO['Sequence with mod'].values[0]
                                                    seq_mod = seq_mod.replace('(Glu->pyro-Glu)','(pyroGlu)')
                                                    seq_mod = seq_mod.replace('(Gln->pyro-Glu)','(pyroGlu)')
                                                    
                                                    fragment_report_path = peptide_report_output + '\\fragment_matches\\' + seq_mod + '_' + str(scan_oi) + '_fragment_report.csv'
                                                    fragment_report = pd.read_csv(fragment_report_path)
                                                    prelim_value_AEIO = fragment_report['Fragment error (Da)'].min()
                                                    
                                                    scan_store.append(scan_oi)
                                                    seq_store.append(seq_oi)
                                                    min_val_store.append(prelim_value_AEIO)
                                                
                                                scan_seq_len = pd.DataFrame()
                                                scan_seq_len['Sequence'] = seq_store
                                                scan_seq_len['Scan'] = scan_store
                                                scan_seq_len['Min Frag Err'] = min_val_store
                                                
                                                value_AEIO = scan_seq_len['Min Frag Err'].min()
                                                criteria_AEIO2 = scan_seq_len[scan_seq_len['Min Frag Err'] == value_AEIO]
                                                criteria_str_AEIO2 = criteria_AEIO2['Sequence'].values.tolist()
                                                
                                                filter_AEIO = filter_AEI[filter_AEI['Sequence'].isin(criteria_str_AEIO2)]
                                                
                                                if len(filter_AEIO) == 1:
                                                    filter_AEIO['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(-) H(-) O(+)'
                                                    final_psm = pd.concat([final_psm,filter_AEIO])
                                                if len(filter_AEIO) > 1:
                                                    filter_AEION = filter_AEIO.drop_duplicates(subset='Sequence')
                                                    if len(filter_AEION) == 1:
                                                        filter_AEION['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(-) H(-) O(++) N(+)'
                                                        final_psm = pd.concat([final_psm,filter_AEION])
                                                    if len(filter_AEION) > 1:
                                                        #raise ValueError('Check this')
                                                        psm_entry = filter_AEION.sample(n=1, random_state=1)
                                                        psm_entry['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(-) H(-) O(++) N(-)'
                                                        final_psm = pd.concat([final_psm,psm_entry])
                                                peptide_list_same_len = criteria_AEILNH['Sequence'].values.tolist()
                                                if len(peptide_list_same_len) == 2:
                                                    seq1 = peptide_list_same_len[0]
                                                    seq2 = peptide_list_same_len[1]
                                                    differing_indexes = indexes = [i for i in range(len(seq1)) if seq1[i] != seq2[i]]
                                                    s = pd.Series(differing_indexes)
                                                    sgc = s.groupby(s.diff().ne(1).cumsum()).transform('count')
                                                    result = s[sgc == sgc.max()].tolist()
                                                    if len(result) <= max_adjacent_swapped_AA:
                                                        filter_AEILNH = criteria_AEILNH[criteria_AEILNH['Sequence'].isin(peptide_list_same_len)]
                                                        filter_AEILNH['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(++) (N++) H(+)'
                                                        final_psm = pd.concat([final_psm,filter_AEILNH])
                                                    if len(result) > max_adjacent_swapped_AA:
                                                        scan_store = []
                                                        seq_store = []
                                                        min_val_store = []
                                                        for line in range(0,len(filter_AEI)):
                                                            criteria_AEILNO = filter_AEI.iloc[[line]]
                                                            scan_oi = criteria_AEILNO['Scan'].values[0]
                                                            seq_oi = criteria_AEILNO['Sequence'].values[0]
                                                            
                                                            seq_mod = criteria_AEILNO['Sequence with mod'].values[0]
                                                            seq_mod = seq_mod.replace('(Glu->pyro-Glu)','(pyroGlu)')
                                                            seq_mod = seq_mod.replace('(Gln->pyro-Glu)','(pyroGlu)')
                                                            
                                                            fragment_report_path = peptide_report_output + '\\fragment_matches\\' + seq_mod + '_' + str(scan_oi) + '_fragment_report.csv'
                                                            fragment_report = pd.read_csv(fragment_report_path)
                                                            prelim_value_AEILNO = fragment_report['Fragment error (Da)'].min()
                                                            
                                                            scan_store.append(scan_oi)
                                                            seq_store.append(seq_oi)
                                                            min_val_store.append(prelim_value_AEILNO)
                                                        
                                                        scan_seq_len = pd.DataFrame()
                                                        scan_seq_len['Sequence'] = seq_store
                                                        scan_seq_len['Scan'] = scan_store
                                                        scan_seq_len['Min Frag Err'] = min_val_store
                                                        
                                                        value_AEILNO = scan_seq_len['Min Frag Err'].min()
                                                        criteria_AEILNO2 = scan_seq_len[scan_seq_len['Min Frag Err'] == value_AEILNO]
                                                        criteria_str_AEILNO2 = criteria_AEILNO2['Sequence'].values.tolist()
                                                        
                                                        filter_AEILNO = filter_AEI[filter_AEI['Sequence'].isin(criteria_str_AEILNO2)]
                                                        
                                                        if len(filter_AEIO) == 1:
                                                            filter_AEILNO['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(++) (N++) H(-) O(+)'
                                                            final_psm = pd.concat([final_psm,filter_AEILNO])
                                                        if len(filter_AEIO) > 1:
                                                            #raise ValueError('Check this')
                                                            psm_entry = filter_AEILNO.sample(n=1, random_state=1)
                                                            psm_entry['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(++) (N++) H(-) O(++?)'
                                                            final_psm = pd.concat([final_psm,psm_entry])
                                                if len(peptide_list_same_len) != 2:
                                                    scan_store = []
                                                    seq_store = []
                                                    min_val_store = []
                                                    for line in range(0,len(filter_AEI)):
                                                        criteria_AEIO = filter_AEI.iloc[[line]]
                                                        scan_oi = criteria_AEIO['Scan'].values[0]
                                                        seq_oi = criteria_AEIO['Sequence'].values[0]
                                                        
                                                        seq_mod = criteria_AEIO['Sequence with mod'].values[0]
                                                        seq_mod = seq_mod.replace('(Glu->pyro-Glu)','(pyroGlu)')
                                                        seq_mod = seq_mod.replace('(Gln->pyro-Glu)','(pyroGlu)')
                                                        
                                                        fragment_report_path = peptide_report_output + '\\fragment_matches\\' + seq_mod + '_' + str(scan_oi) + '_fragment_report.csv'
                                                        fragment_report = pd.read_csv(fragment_report_path)
                                                        prelim_value_AEIO = fragment_report['Fragment error (Da)'].min()
                                                        
                                                        scan_store.append(scan_oi)
                                                        seq_store.append(seq_oi)
                                                        min_val_store.append(prelim_value_AEIO)
                                                    
                                                    scan_seq_len = pd.DataFrame()
                                                    scan_seq_len['Sequence'] = seq_store
                                                    scan_seq_len['Scan'] = scan_store
                                                    scan_seq_len['Min Frag Err'] = min_val_store
                                                    
                                                    value_AEIO = scan_seq_len['Min Frag Err'].min()
                                                    criteria_AEIO2 = scan_seq_len[scan_seq_len['Min Frag Err'] == value_AEIO]
                                                    criteria_str_AEIO2 = criteria_AEIO2['Sequence'].values.tolist()
                                                    
                                                    filter_AEIO = filter_AEI[filter_AEI['Sequence'].isin(criteria_str_AEIO2)]
                                                    
                                                    if len(filter_AEIO) == 1:
                                                        filter_AEIO['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(-) H(-) O(+)'
                                                        final_psm = pd.concat([final_psm,filter_AEIO])
                                                    if len(filter_AEIO) > 1:
                                                        filter_AEION = filter_AEIO.drop_duplicates(subset='Sequence')
                                                        if len(filter_AEION) == 1:
                                                            filter_AEION['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(-) H(-) O(++) N(+)'
                                                            final_psm = pd.concat([final_psm,filter_AEION])
                                                        if len(filter_AEION) > 1:
                                                            #raise ValueError('Check this')
                                                            psm_entry = filter_AEION.sample(n=1, random_state=1)
                                                            psm_entry['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(-) H(-) O(++) N(-)'
                                                            final_psm = pd.concat([final_psm,psm_entry])
                                        if len(filter_AEILN) != len(criteria_AEILNH): 
                                             #raise ValueError('Check this') #no impact, validated 6/28/23
                                             psm_entry = filter_AEI.sample(n=1, random_state=1)
                                             psm_entry['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(-) H(-?)'
                                             final_psm = pd.concat([final_psm,psm_entry])
    
        
                            else:
                                filter_AEI['Seq_Len']  = filter_AEI['Sequence'].str.len()
                                max_seq_len = filter_AEI['Seq_Len'].max()
                                criteria_AEIH = filter_AEI[filter_AEI['Seq_Len'] == max_seq_len]
                                if len(filter_AEI) == len(criteria_AEIH): #filter for peptides of same length
                                    peptide_list_same_len = criteria_AEIH['Sequence'].values.tolist()
                                    if len(peptide_list_same_len) == 2:
                                        seq1 = peptide_list_same_len[0]
                                        seq2 = peptide_list_same_len[1]
                                        differing_indexes = indexes = [i for i in range(len(seq1)) if seq1[i] != seq2[i]]
                                        s = pd.Series(differing_indexes)
                                        sgc = s.groupby(s.diff().ne(1).cumsum()).transform('count')
                                        result = s[sgc == sgc.max()].tolist()
                                        if len(result) <= max_adjacent_swapped_AA:
                                            filter_AEIH = criteria_AEIH[criteria_AEIH['Sequence'].isin(peptide_list_same_len)]
                                            filter_AEIH['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(-) H(+)'
                                            final_psm = pd.concat([final_psm,filter_AEIH])
                                        if len(result) > max_adjacent_swapped_AA:
                                            scan_store = []
                                            seq_store = []
                                            min_val_store = []
                                            for line in range(0,len(filter_AEI)):
                                                criteria_AEIO = filter_AEI.iloc[[line]]
                                                scan_oi = criteria_AEIO['Scan'].values[0]
                                                seq_oi = criteria_AEIO['Sequence'].values[0]
                                                
                                                seq_mod = criteria_AEIO['Sequence with mod'].values[0]
                                                seq_mod = seq_mod.replace('(Glu->pyro-Glu)','(pyroGlu)')
                                                seq_mod = seq_mod.replace('(Gln->pyro-Glu)','(pyroGlu)')
                                                
                                                fragment_report_path = peptide_report_output + '\\fragment_matches\\' + seq_mod + '_' + str(scan_oi) + '_fragment_report.csv'
                                                fragment_report = pd.read_csv(fragment_report_path)
                                                prelim_value_AEIO = fragment_report['Fragment error (Da)'].min()
                                                
                                                scan_store.append(scan_oi)
                                                seq_store.append(seq_oi)
                                                min_val_store.append(prelim_value_AEIO)
                                            
                                            scan_seq_len = pd.DataFrame()
                                            scan_seq_len['Sequence'] = seq_store
                                            scan_seq_len['Scan'] = scan_store
                                            scan_seq_len['Min Frag Err'] = min_val_store
                                            
                                            value_AEIO = scan_seq_len['Min Frag Err'].min()
                                            criteria_AEIO2 = scan_seq_len[scan_seq_len['Min Frag Err'] == value_AEIO]
                                            criteria_str_AEIO2 = criteria_AEIO2['Sequence'].values.tolist()
                                            
                                            filter_AEIO = filter_AEI[filter_AEI['Sequence'].isin(criteria_str_AEIO2)]
                                            
                                            if len(filter_AEIO) == 1:
                                                filter_AEIO['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(-) H(-) O(+)'
                                                final_psm = pd.concat([final_psm,filter_AEIO])
                                            if len(filter_AEIO) > 1:
                                                #raise ValueError('Check this')
                                                psm_entry = filter_AEIO.sample(n=1, random_state=1)
                                                psm_entry['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(-) H(-) O(++?)'
                                                final_psm = pd.concat([final_psm,psm_entry])
                                    if len(peptide_list_same_len) != 2:
                                        scan_store = []
                                        seq_store = []
                                        min_val_store = []
                                        for line in range(0,len(filter_AEI)):
                                            criteria_AEIO = filter_AEI.iloc[[line]]
                                            scan_oi = criteria_AEIO['Scan'].values[0]
                                            seq_oi = criteria_AEIO['Sequence'].values[0]
                                            
                                            seq_mod = criteria_AEIO['Sequence with mod'].values[0]
                                            seq_mod = seq_mod.replace('(Glu->pyro-Glu)','(pyroGlu)')
                                            seq_mod = seq_mod.replace('(Gln->pyro-Glu)','(pyroGlu)')
                                            
                                            fragment_report_path = peptide_report_output + '\\fragment_matches\\' + seq_mod + '_' + str(scan_oi) + '_fragment_report.csv'
                                            fragment_report = pd.read_csv(fragment_report_path)
                                            prelim_value_AEIO = fragment_report['Fragment error (Da)'].min()
                                            
                                            scan_store.append(scan_oi)
                                            seq_store.append(seq_oi)
                                            min_val_store.append(prelim_value_AEIO)
                                        
                                        scan_seq_len = pd.DataFrame()
                                        scan_seq_len['Sequence'] = seq_store
                                        scan_seq_len['Scan'] = scan_store
                                        scan_seq_len['Min Frag Err'] = min_val_store
                                        
                                        value_AEIO = scan_seq_len['Min Frag Err'].min()
                                        criteria_AEIO2 = scan_seq_len[scan_seq_len['Min Frag Err'] == value_AEIO]
                                        criteria_str_AEIO2 = criteria_AEIO2['Sequence'].values.tolist()
                                        
                                        filter_AEIO = filter_AEI[filter_AEI['Sequence'].isin(criteria_str_AEIO2)]
                                        
                                        if len(filter_AEIO) == 1:
                                            filter_AEIO['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(-) H(-) O(+)'
                                            final_psm = pd.concat([final_psm,filter_AEIO])
                                        if len(filter_AEIO) > 1:
                                            filter_AEION = filter_AEIO.drop_duplicates(subset='Sequence')
                                            if len(filter_AEION) == 1:
                                                filter_AEION['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(-) H(-) O(++) N(+)'
                                                final_psm = pd.concat([final_psm,filter_AEION])
                                            if len(filter_AEION) > 1:
                                                if len(filter_AEION) == 2:
                                                    criteria_str_AEIONP = filter_AEION['Sequence'].values.tolist()
                                                    string1 = criteria_str_AEIONP[0]
                                                    string2 = criteria_str_AEIONP[1]
                                                    num_dif = dif_compare(string1,string2)
                                                    if num_dif <= num_sub_AAs:
                                                        filter_AEIONP = filter_AEION
                                                        filter_AEIONP['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(-) H(-) O(++) N(+) P(+)'
                                                        final_psm = pd.concat([final_psm,filter_AEIONP])
                                                    if num_dif > num_sub_AAs:
                                                        #raise ValueError('Check this')
                                                        psm_entry = filter_AEION.sample(n=1, random_state=1)
                                                        psm_entry['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(-) H(-) O(++) N(-) P(-)'
                                                        final_psm = pd.concat([final_psm,psm_entry])
                                                if len(filter_AEION)> 2:
                                                    #raise ValueError('Check this')
                                                    psm_entry = filter_AEION.sample(n=1, random_state=1)
                                                    psm_entry['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(-) H(-) O(++) N(-) P(-)'
                                                    final_psm = pd.concat([final_psm,psm_entry])
                                        
                                if len(filter_AEI) != len(criteria_AEIH): 
                                     
                                     value_AEIG = filter_AEI['Correlation value'].max()
                                     filter_AEIG = filter_AEI[filter_AEI['Correlation value'] == value_AEIG]
                                     if len(filter_AEIG) == 1:
                                         filter_AEIG['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(-) H(-) G(+)'
                                         final_psm = pd.concat([final_psm,filter_AEIG])
                                     if len(filter_AEIG) > 1:     
                                         #raise ValueError('Check this')
                                         psm_entry = filter_AEIG.sample(n=1, random_state=1)
                                         psm_entry['Step assigned'] = 'filter A(-) E(-) I(++) F(-) L(-) H(-) G(-)'
                                         final_psm = pd.concat([final_psm,psm_entry])
    
        
        first_round_psm_no_dups = final_psm.drop_duplicates(subset='Sequence')
        first_round_psm_no_dups_target = first_round_psm_no_dups[first_round_psm_no_dups['Sequence'].isin(target_list)]
        first_round_psm_no_dups_decoy = first_round_psm_no_dups[~first_round_psm_no_dups['Sequence'].isin(target_list)]
        
        psm_total_target = final_psm[final_psm['Sequence'].isin(target_list)]
        psm_total_decoy = final_psm[~final_psm['Sequence'].isin(target_list)]
        
        print(directory)
        print('# PSMs: ',len(final_psm))
        print('# Unique IDs: ',len(first_round_psm_no_dups))
        print('# Unique Target IDs: ',len(first_round_psm_no_dups_target))
        print('# Unique Decoy IDs: ',len(first_round_psm_no_dups_decoy))
        print('Ratio Target:Decoy IDs: ', (len(first_round_psm_no_dups_target)/len(first_round_psm_no_dups_decoy)))
        print('# Target PSMs: ', (len(psm_total_target)))
        print('# Decoy PSMs: ', (len(psm_total_decoy)))
        print('FDR @ PSM level: ', (len(psm_total_decoy)/len(psm_total_target)))
        
        final_psm_target = final_psm[final_psm['Sequence'].isin(target_list)]
        final_psm_decoy = final_psm[~final_psm['Sequence'].isin(target_list)]
        
        final_psm_target['Status'] = 'Target'
        final_psm_decoy['Status'] = 'Decoy'
        
        final_psm_output = pd.concat([final_psm_target,final_psm_decoy])
        
        output_path = peptide_report_output + '\\final_psm_report_out.csv'
        with open(output_path,'w',newline='') as filec:
                writerc = csv.writer(filec)
                final_psm_output.to_csv(filec,index=False)
                
        
        return final_psm_output
        
        text_file_path = db_search_parent_directory + '\\PSM_summary.txt'
        file1 = open(text_file_path, "a")
        L = ['\n\n\n' + (directory) +
        ('\n# PSMs: ' + str(len(final_psm))) +
        ('\n# Unique IDs: ' + str(len(first_round_psm_no_dups))) +
        ('\n# Unique Target IDs: ' + str(len(first_round_psm_no_dups_target))) +
        ('\n# Unique Decoy IDs: ' + str(len(first_round_psm_no_dups_decoy))) + 
        ('\nRatio Target:Decoy IDs: ' + str((len(first_round_psm_no_dups_target)/len(first_round_psm_no_dups_decoy)))) +
        ('\n# Target PSMs: ' + str((len(psm_total_target)))) +
        ('\n# Decoy PSMs: ' + str((len(psm_total_decoy)))) +
        ('\nFDR @ PSM level: ' + str((len(psm_total_decoy)/len(psm_total_target))))]
        file1.writelines(L)
        file1.close()
    
    subject = 'Your code has finished running'
    text = 'PSM assignment for all samples has finished running'
    content = 'Subject: %s\n\n%s' % (subject, text)
    mail = smtplib.SMTP('smtp.gmail.com',587)
    mail.ehlo()
    mail.starttls()
    mail.login('lingjun.li.notifications@gmail.com','eabtnjwaikdssdtd')
    mail.sendmail('lingjun.li.notifications@gmail.com','lawashburn@wisc.edu',content) 
    mail.close()