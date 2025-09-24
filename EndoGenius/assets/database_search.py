# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 16:03:08 2023

@author: lawashburn
"""

import csv
import pandas as pd
import re
import os
from itertools import permutations
from Bio.SeqIO.FastaIO import SimpleFastaParser
import random
import collections
import time
import numpy as np
from scipy import signal
from datetime import datetime
import scipy
from scipy.spatial import distance
import numpy as np
from pyopenms import *
import smtplib

def launch_db_search_pt1(predefined_db_path,output_parent_directory,raw_file_formatted_path,precursor_error_cutoff,
                         fragment_error_cutoff,min_mz,min_intensity,standard_err_percent,amidation,oxidation_M_status,
                         pyroglu_E_status,pyroglu_Q_status,sulfo_Y_status,max_modifications):
   
    backslash_index = raw_file_formatted_path.rfind('\\')
    base_file_path = raw_file_formatted_path[0:(backslash_index)]
    
    # base_file_path = r"C:\Users\lawashburn\Documents\DBpep_v2\validation\Nhu_RAW_files"
    # output_parent_directory = r"C:\Users\lawashburn\Documents\DB_pep_validation\DB_search_w_mod_20230629"
    
    mzml_path = raw_file_formatted_path.replace('.txt','.mzML')

    raw_conv_mzml_storage = [[mzml_path,raw_file_formatted_path]]

    # raw_conv_mzml_storage = [[r"C:\Users\lawashburn\Documents\DBpep_v2\validation\Nhu_RAW_files\2021_0817_SG_3.mzML",
    #                           r"C:\Users\lawashburn\Documents\DBpep_v2\validation\Nhu_RAW_files\2021_0817_SG_3_formatted.txt"],
    #                          [r"C:\Users\lawashburn\Documents\DBpep_v2\validation\Nhu_RAW_files\2021_0817_TG_1.mzML",
    #                           r"C:\Users\lawashburn\Documents\DBpep_v2\validation\Nhu_RAW_files\2021_0817_TG_1_formatted.txt"],
    #                          [r"C:\Users\lawashburn\Documents\DBpep_v2\validation\Nhu_RAW_files\2021_0817_TG_2.mzML",
    #                           r"C:\Users\lawashburn\Documents\DBpep_v2\validation\Nhu_RAW_files\2021_0817_TG_2_formatted.txt"],
    #                          [r"C:\Users\lawashburn\Documents\DBpep_v2\validation\Nhu_RAW_files\2021_0817_TG_3.mzML",
    #                           r"C:\Users\lawashburn\Documents\DBpep_v2\validation\Nhu_RAW_files\2021_0817_TG_3_formatted.txt"]]
    
    h_mass = 1.00784
    # precursor_error_cutoff = 20 #ppm
    # fragment_error_cutoff = 0.02
    # min_mz = 50
    # min_intensity = 1000
    # standard_err_percent = 0.1
    spectra_segments = 50
    
    # amidation = True
    # oxidation_M_status = True
    # pyroglu_E_status = True
    # pyroglu_Q_status = True
    # sulfo_Y_status = True
    # max_modifications = 3
    
    ### generating database of mods selected ###
    mods = []
    mod_dict = {}
    if oxidation_M_status == 1:
        mods.append('M(Oxidation)')
        mod_dict['M'] = 'M(Oxidation)'
    if pyroglu_E_status == 1:
        mods.append('E(Glu->pyro-Glu)')
        mod_dict['E'] = 'E(Glu->pyro-Glu)'
    if pyroglu_Q_status == 1:
        mods.append('Q(Gln->pyro-Glu)')
        mod_dict['Q'] = 'Q(Gln->pyro-Glu)'
    if sulfo_Y_status == 1:
        mods.append('Y(Sulfo)')
        mod_dict['Y'] = 'Y(Sulfo)'
    else:
        pass
    
    modded_aas = []
    for a in mods:
        modded_aas.append(a[0])
    
    ### Theoretical fragment calculator ###
    proton_mass = 1.00727646688
    charge = 1 #fragment charge
    H = 1.0078250352
    O = 15.99491463
    C = 12.0000000
    N = 14.003074
    P = 30.973762
    S = 31.9720707
    
    aa_masses = {
        'G' : C*2  + H*3  + N   + O,
        'A' : C*3  + H*5  + N   + O,
        'S' : C*3  + H*5  + N   + O*2,
        'P' : C*5  + H*7  + N   + O,
        'V' : C*5  + H*9  + N   + O,
        'T' : C*4  + H*7  + N   + O*2,
        'C' : C*3  + H*5  + N   + O   + S,
        'L' : C*6  + H*11 + N   + O,
        'I' : C*6  + H*11 + N   + O,
        'N' : C*4  + H*6  + N*2 + O*2,
        'D' : C*4  + H*5  + N   + O*3,
        'Q' : C*5  + H*8  + N*2 + O*2,
        'K' : C*6  + H*12 + N*2 + O ,
        'E' : C*5  + H*7  + N   + O*3 ,
        'M' : C*5  + H*9  + N   + O   + S ,
        'H' : C*6  + H*7  + N*3 + O ,
        'F' : C*9  + H*9  + N   + O ,
        'R' : C*6  + H*12 + N*4 + O ,
        'Y' : C*9  + H*9  + N   + O*2 ,
        'W' : C*11 + H*10 + N*2 + O ,
        'O' : C*5  + H*12 + N*2 + O*2,
        'C(Pyro-glu)' : C*3  + H * 2 + O + S,
        'Q(Gln->pyro-Glu)' : C*5  + H*5  + N + O*2,
        'E(Glu->pyro-Glu)' : C*5  + H*4 + O*3,
        'M(Oxidation)' : C*5  + H*9  + N   + O*2   + S,
        'Y(Sulfo)' :  C*9  + H*9  + N   + O*5 + S
        }
    
    termini = {'Standard' : H * 2 + O,
    '(Amidated)' : N + H * 3}
    PTMs = {'C(Pyro-glu)' : H * -3 - N,
            'Q(Gln->pyro-Glu)' : H * -3 - N,
            'E(Glu->pyro-Glu)' : H * -3 - N,
            'M(Oxidation)' : O,
            'Y(Sulfo)' : S + O * 3
            }
    adducts = {
        'H2O' : H * 2 + O,
        'NH3' : N + H * 3}
    
    #%%
    def seq_coverage_calc(merge_fragment_match_filtered,ion_report,scan_to_report,peptide):
        frag_merge_match_filtered = merge_fragment_match_filtered[merge_fragment_match_filtered['Fragment error (Da)'] <= fragment_error_cutoff]
        if len(frag_merge_match_filtered)>0:
            frag_merge_match_filtered['ion_format'] = frag_merge_match_filtered['ion'].str.extract('(\d+)', expand=False)
            ion_report['ion_format'] = ion_report['ion'].str.extract('(\d+)', expand=False)
            merge_fragment_match_filtered_ion = frag_merge_match_filtered.drop_duplicates(subset=['ion_format'],keep='first').reset_index(drop=True)
            ion_report_filtered_ion = ion_report.drop_duplicates(subset=['ion_format'],keep='first').reset_index(drop=True)
            seq_coverage = (len(merge_fragment_match_filtered_ion)/len(ion_report_filtered_ion))*100
            
            seq_data = {'Sequence':[peptide],
                        'Scan':[scan_to_report],
                        'Sequence coverage':[seq_coverage]}
            
            seq_coverage_rep = pd.DataFrame(seq_data)
            return seq_coverage_rep
        else:
            seq_data = {'Sequence':[peptide],
                        'Scan':[scan_to_report],
                        'Sequence coverage':[0]}
            
            seq_coverage_rep = pd.DataFrame(seq_data)
            return seq_coverage_rep
    #     #%%
    
    def hyperscore_calc(ss,peptide):        
        ss = int(ss)
        tsg = TheoreticalSpectrumGenerator()
        thspec = MSSpectrum()
        p = Param()
        p.setValue("add_metainfo", "true")
        tsg.setParameters(p)
        peptide = AASequence.fromString(peptide)
        tsg.getSpectrum(thspec, peptide, 1, 1)
    
        spectrum_of_interest = e[ss-1]
        print("Spectrum native id", spectrum_of_interest.getNativeID())
        mz, i = spectrum_of_interest.get_peaks()
        peaks = [(mz, i) for mz, i in zip(mz, i) if i > min_intensity and mz > min_mz]
    
        hscore = HyperScore()
        fragment_mass_tolerance = precursor_error_cutoff
        is_tol_in_ppm = True
        result = hscore.compute(
            fragment_mass_tolerance, is_tol_in_ppm, spectrum_of_interest, thspec)
    
        return result
    #%%
    def check_termini_Key(dict, key):
        if key in dict.keys():
            return dict[key]
        else:
            return termini['Standard']
    
    def check_PTM_Key(dict, key):
        if key in dict.keys():
            return dict[key]
        else:
            return 0
    
    def check_termini(pot_mod_peptide):
        if '(' in pot_mod_peptide:
            term_start = (pot_mod_peptide.rindex('('))
            termini_ID = pot_mod_peptide[term_start:]
            termini_mass_change = check_termini_Key(termini, termini_ID)
            return termini_mass_change
        else:
            return termini['Standard']
    
    def check_PTM(pot_mod_peptide):
        number_of_mods = pot_mod_peptide.count('(')
        if number_of_mods > 0:
            current_peptide = []
            mass_change_collection = []
            current_peptide.append(pot_mod_peptide)
            for a in range(0,number_of_mods):
                peptide_less_mod = current_peptide[-1]
                ptm_start = (peptide_less_mod.index('('))-1
                ptm_end = (peptide_less_mod.index(')'))+1
                ptm_ID = peptide_less_mod[ptm_start:ptm_end]
                ptm_mass_change = check_PTM_Key(PTMs, ptm_ID)
                mass_change_collection.append(ptm_mass_change)
                peptide_less_mod2 = peptide_less_mod[:ptm_start] + peptide_less_mod[ptm_end:]
                current_peptide.append(peptide_less_mod2)
                
            ptm_mass_change = sum(mass_change_collection)
            return ptm_mass_change
        else:
            ptm_mass_change = 0
            return ptm_mass_change
    
    def list_of_residues(pot_mod_peptide):
        list_of_res = []
        pep_update = []
        pep_update.append(pot_mod_peptide)
        no_mods = pot_mod_peptide.count('(')
        if no_mods > 1:
            for c in range(0,no_mods+1):
                pep_of_interest = pep_update[-1]
                if '(' in pep_of_interest:
                    first_ptm_start = pep_of_interest.index('(')
                    first_ptm_end = pep_of_interest.index(')')
    
                    first_residues = pep_of_interest[:(first_ptm_start-1)]
                    for a in first_residues:
                        list_of_res.append(a)
                    ptm_residue = pep_of_interest[(first_ptm_start-1):(first_ptm_end+1)]
                    list_of_res.append(ptm_residue)
                    remaining_pep = pep_of_interest[(first_ptm_end+1):]
                    pep_update.append(remaining_pep)  
                else:
                    for d in pep_of_interest:
                        list_of_res.append(d)
        elif no_mods == 1:
            for c in range(0,1):
                pep_of_interest = pep_update[-1]
                if '(' in pep_of_interest:
                    first_ptm_start = pep_of_interest.index('(')
                    first_ptm_end = pep_of_interest.index(')')
                    if first_ptm_start == 1:
                        ptm_residue =  pep_of_interest[0] + (pep_of_interest[(first_ptm_start):(first_ptm_end+1)])
                        list_of_res.append(ptm_residue)
                        remaining_pep = pep_of_interest[(first_ptm_end+1):]
                        for d in remaining_pep:
                            list_of_res.append(d)
                    if first_ptm_start != 1:
                        first_residues = pep_of_interest[:(first_ptm_start-1)]
                        for a in first_residues:
                            list_of_res.append(a)
                        ptm_residue = pep_of_interest[(first_ptm_start-1):(first_ptm_end+1)]
                        list_of_res.append(ptm_residue)
                        remaining_pep = pep_of_interest[(first_ptm_end+1):]
                        for d in remaining_pep:
                            list_of_res.append(d)              
                else:
                    for d in pep_of_interest:
                        list_of_res.append(d) 
        elif no_mods == 0:
            for c in pot_mod_peptide:
                list_of_res.append(c)
        return list_of_res
    ### end of theoretical fragment calculator ###
    ### start of monoisotopic mass calculator ###
    def monoisotopic_mass_calculator(peptide_from_fasta):
            plain_peptide = re.sub("[\(\[].*?[\)\]]", "",peptide_from_fasta)
            res_list_for_fragment = list_of_residues(peptide_from_fasta)
            mass_of_residues = []
            for residue in plain_peptide:
                residue_mass = aa_masses[residue]
                mass_of_residues.append(residue_mass)
            peptide_mass = (sum(mass_of_residues)) + check_termini(peptide_from_fasta) + check_PTM(peptide_from_fasta)
            mass_to_charge = (peptide_mass + (proton_mass * charge))/charge
            return mass_to_charge
    ### end of monoisotopic mass calculator ###
    
    ## Start of theoretical spectra generator ##
    def theoretical_spectra_generator(peptide_to_check):
        ###pull theoretical masses###
        plain_peptide = re.sub("[\(\[].*?[\)\]]", "",peptide_to_check) #removes any modification for mass calculations
        res_list_for_fragment = list_of_residues(peptide_to_check)
        mass_of_residues = []
        for residue in plain_peptide:
            residue_mass = aa_masses[residue]
            mass_of_residues.append(residue_mass)
        peptide_mass = (sum(mass_of_residues)) + check_termini(peptide_to_check) + check_PTM(peptide_to_check) #calculated MH mass
        mass_to_charge = (peptide_mass + (proton_mass * charge))/charge #translates the MH mass to m/z for each charge of interest
    
        num_ions = len(plain_peptide)-1 #number of expected fragment ions is one less than the number of AAs in the sequence
    
        b_ions = []
        b_ion_name = []
        
        y_ions = []
        y_ion_name = []
        
        for a in range(0,num_ions):
            
            residue_identity = res_list_for_fragment[a]
            if len(b_ions) == 0:
                ion_mass = aa_masses[residue_identity]
                ion_mz = ion_mass + proton_mass
                b_ions.append(ion_mz)
                ion_name = 'b' + str(a+1)
                b_ion_name.append(ion_name)
            
            elif len(b_ions) > 0:
                ion_mass = (aa_masses[residue_identity]) + b_ions[-1]
                b_ions.append(ion_mass)
                ion_name = 'b' + str(a+1)
                b_ion_name.append(ion_name)
        
        for b in (range(0,num_ions)):
            residue_identity = res_list_for_fragment[b]
            if len(y_ions) == 0:
                ion_mass = mass_to_charge - aa_masses[residue_identity]
                y_ions.append(ion_mass)
                ion_name = 'y' + str((num_ions-b))
                y_ion_name.append(ion_name)
            elif len(y_ions) > 0:
                ion_mass = y_ions[-1] - aa_masses[residue_identity]
                y_ions.append(ion_mass)
                ion_name = 'y' + str((num_ions-b))
                y_ion_name.append(ion_name)
    
        b_ions_report = pd.DataFrame()
        b_ions_report['ion'] = b_ion_name
        b_ions_report['mass'] = b_ions
        
        b_ions_water_adduct = pd.DataFrame()
        b_ions_water_adduct['ion'] = b_ions_report['ion'] + '-H2O'
        b_ions_water_adduct['mass'] = b_ions_report['mass'] - adducts['H2O']
        
        b_ions_ammonia_adduct = pd.DataFrame()
        b_ions_ammonia_adduct['ion'] = b_ions_report['ion'] + '-NH3'
        b_ions_ammonia_adduct['mass'] = b_ions_report['mass'] - adducts['NH3']
        
        y_ions_report = pd.DataFrame()
        y_ions_report['ion'] = y_ion_name
        y_ions_report['mass'] = y_ions
        
        y_ions_ammonia_adduct = pd.DataFrame()
        y_ions_ammonia_adduct['ion'] = y_ions_report['ion'] + '-NH3'
        y_ions_ammonia_adduct['mass'] = y_ions_report['mass'] - adducts['NH3']
    
        y_ions_water_adduct = pd.DataFrame()
        y_ions_water_adduct['ion'] = y_ions_report['ion'] + '-H2O'
        y_ions_water_adduct['mass'] = y_ions_report['mass'] - adducts['H2O']
        
        ion_report = pd.DataFrame()
        ion_report = pd.concat([ion_report,b_ions_report])
        ion_report = pd.concat([ion_report,y_ions_report])
        ion_report = pd.concat([ion_report,b_ions_water_adduct])
        ion_report = pd.concat([ion_report,b_ions_ammonia_adduct])
        ion_report = pd.concat([ion_report,y_ions_ammonia_adduct])
        ion_report = pd.concat([ion_report,y_ions_water_adduct])
        ion_report = ion_report.drop_duplicates()
        
        ion_report = ion_report.rename(columns={'mass':'Fragment theoretical monoisotopic mass'})
        return ion_report
    
    ###Definitions###
    def raw_file_detail_extraction(raw_file_path):
        raw_file_sample_name1 = raw_converter_path.replace(base_file_path,'')
        raw_file_sample_name2 = raw_file_sample_name1.replace('_formatted','')
        raw_file_sample_name3 = raw_file_sample_name2.replace('\\','')
        sample_name = raw_file_sample_name3.replace('.txt','')
        output_folder = output_parent_directory+'\\'+sample_name
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        ### generate output folder ###
        return sample_name, output_folder
    
    def raw_file_data_extraction(raw_file_path):
        raw_converter = pd.read_csv(raw_converter_path, sep=",",skiprows=[0], names= ["m/z","resolution","charge","intensity", "MS2", "scan_number","precursor_charge",'null','Sample','Identifier','Iteration'])
        raw_converter = raw_converter.rename(columns={'m/z':'Fragment actual m/z',
                                                      'charge': 'Fragment actual charge',
                                                      'intensity':'Fragment actual intensity',
                                                      'MS2':'Precursor actual m/z',
                                                      "precursor_charge":'Precursor actual charge',
                                                      'scan_number':'Scan'})
    
        raw_converter['Fragment actual charge'] = raw_converter['Fragment actual charge'].replace(to_replace=0,value=1) #assume z=0 is z=1   
        exp_precursor = raw_converter.drop_duplicates() 
        exp_precursor = exp_precursor.copy()
        exp_precursor['Precursor actual monoisotopic mass'] =  ((exp_precursor['Precursor actual m/z']) * 
                                                                (exp_precursor['Precursor actual charge']))-(h_mass*(exp_precursor['Precursor actual charge']))
        return exp_precursor
    
    
    fasta_w_mass = pd.read_csv(predefined_db_path)
    
    for file_category in raw_conv_mzml_storage:
        
        raw_converter_path = file_category[1]
        mzml_path_input = file_category[0]
        rounds_number = [1]
        unique_IDS_number = []
        details = raw_file_detail_extraction(raw_converter_path)
        sample_name = details[0]
    
        peptide_report_output = details[1]
        exp_precursor = raw_file_data_extraction(raw_converter_path) 
        exp_precursor = exp_precursor.drop(columns=['resolution','Sample','Identifier','Iteration','null'])
        # fasta_w_mass = fasta_w_mass[fasta_w_mass["Sequence"].str.contains("\(") == False]
        
        #precursor AMM#
        precursor_temp_cutoff = precursor_error_cutoff*3 #rough estimate of ppm to Da to minimize extra search space
     
        if len(fasta_w_mass)<1: #throws an error if database file is empty
            raise ValueError('Database file is empty')
        
        exp_mass_only = exp_precursor.drop(columns=['Fragment actual m/z','Fragment actual charge','Fragment actual intensity','Precursor actual m/z','Precursor actual charge'])
    
        db_mass_only = fasta_w_mass
    
        exp_precursor_sorted = exp_mass_only.drop_duplicates(subset = ['Precursor actual monoisotopic mass','Scan'])
        exp_precursor_sorted = exp_precursor_sorted.sort_values(by = 'Precursor actual monoisotopic mass')
        
        
        db_sorted_filtered = db_mass_only.drop_duplicates(subset='Precursor theoretical monoisotopic mass')
        db_sorted_filtered = db_sorted_filtered.sort_values(by = 'Precursor theoretical monoisotopic mass') #sort database monoisotopic mass and experimental for mergeasof
        
        exp_precursor_sorted_filtered = exp_precursor_sorted.drop_duplicates(subset='Precursor actual monoisotopic mass')
    
        merge_match=pd.DataFrame()
        merge_match_dfs = []
        for a in range(0,len(db_sorted_filtered)):
    
            db_filtered = db_sorted_filtered.iloc[[a]]
    
            merge_match2 = pd.merge_asof(exp_precursor_sorted_filtered,db_sorted_filtered, left_on='Precursor actual monoisotopic mass', 
                                        right_on='Precursor theoretical monoisotopic mass',
                                        tolerance = precursor_temp_cutoff, allow_exact_matches=True,direction='forward') 
    
            merge_match3 = pd.merge_asof(exp_precursor_sorted_filtered,db_sorted_filtered, left_on='Precursor actual monoisotopic mass', 
                                        right_on='Precursor theoretical monoisotopic mass',
                                        tolerance = precursor_temp_cutoff, allow_exact_matches=True,direction='backward') 
    
            merge_match_dfs.append(merge_match2)
            merge_match_dfs.append(merge_match3)
    
        merge_match = pd.concat(merge_match_dfs,ignore_index=True)
    
        
        merge_match = merge_match[merge_match['Sequence'].notna()]
        merge_match = merge_match.drop_duplicates()
        merge_match = merge_match.sort_values(by='Scan')
       
        
        merge_match_merge_again = pd.merge(merge_match, db_mass_only, on='Precursor theoretical monoisotopic mass', how='inner')
        merge_match_merge_again_again = pd.merge(merge_match_merge_again, exp_precursor_sorted, on='Precursor actual monoisotopic mass', how='inner')
        merge_match_count = merge_match_merge_again_again.drop_duplicates(subset=['Scan_y','Sequence_y'])
    
        #end precursor AMM#
        merge_match_count = merge_match_count.drop(columns=['Sequence_x','Scan_x'])
        merge_match_count = merge_match_count.rename(columns={'Sequence_y':'Sequence','Scan_y':'Scan'})
        
        merge_match_count['Precursor error (ppm)'] = ((abs((merge_match_count['Precursor theoretical monoisotopic mass'])-
                                                              (merge_match_count['Precursor actual monoisotopic mass'])))/
                                                          (merge_match_count['Precursor theoretical monoisotopic mass'])) * 1E6
    
        precursor_amm_results = merge_match_count[merge_match_count['Precursor error (ppm)'] <= precursor_error_cutoff]
        
        output_path = peptide_report_output + '\\precursor_AMM_results.csv'
        with open(output_path,'w',newline='') as filec:
                writerc = csv.writer(filec)
                precursor_amm_results.to_csv(filec,index=False)
                
        all_seqs = precursor_amm_results.drop(columns=['Scan','Precursor theoretical monoisotopic mass','Precursor actual monoisotopic mass'])
        # output_path_rep = peptide_report_output + '\\all_seqs1.csv'
        # with open(output_path_rep,'w',newline='') as filec:
        #         writerc = csv.writer(filec)
        #         all_seqs.to_csv(filec,index=False)
        all_seqs = all_seqs.drop_duplicates(subset='Sequence',ignore_index=True)
        #%%
        unfiltered_psms = []
        # output_path_rep = peptide_report_output + '\\all_seqs.csv'
        # with open(output_path_rep,'w',newline='') as filec:
        #         writerc = csv.writer(filec)
        #         all_seqs.to_csv(filec,index=False)
    
        for a in range(0,len(precursor_amm_results)):
            candidates_filtered_AMM = precursor_amm_results.iloc[[a]]
            # print(candidates_filtered)
            candidate_seq = candidates_filtered_AMM['Sequence'].values[0]
            candidate_scan = candidates_filtered_AMM['Scan'].values[0]
            
            candidates_filtered = all_seqs[all_seqs['Sequence'] == candidate_seq]
            
            ion_report = theoretical_spectra_generator(candidate_seq)
            
            scan_extract = pd.merge(candidates_filtered,precursor_amm_results,on='Sequence')
            fragments_extract = pd.merge(scan_extract,exp_precursor,on=['Scan'])
            
            fragments_extract['Fragment actual monoisotopic mass'] = (fragments_extract['Fragment actual m/z'] * 
                                                                           fragments_extract['Fragment actual charge']) - (h_mass*fragments_extract['Fragment actual charge'])  
            fragments_extract = fragments_extract.sort_values(by='Fragment actual monoisotopic mass')
            
            # scan_extract = pd.merge(candidates_filtered,precursor_amm_results,on='Sequence')
            
            fragments_extract = fragments_extract[fragments_extract['Scan'] == candidate_scan]
            
            # fragments_extract = pd.merge(scan_extract,exp_precursor,on=['Scan'])
            
            # output_path_rep = peptide_report_output + '\\fragments_extract.csv'
            # with open(output_path_rep,'w',newline='') as filec:
            #         writerc = csv.writer(filec)
            #         fragments_extract.to_csv(filec,index=False)
            
            #Begin fragment AMM
            
            
            # output_path_rep = peptide_report_output + '\\fragments_extract2.csv'
            # with open(output_path_rep,'w',newline='') as filec:
            #         writerc = csv.writer(filec)
            #         fragments_extract.to_csv(filec,index=False)
            
            scans_present_raw = fragments_extract['Scan'].values.tolist()
            scans_present = []
            for z in scans_present_raw:
                if z not in scans_present:
                    scans_present.append(z)
            peptide_rep_output_folder = peptide_report_output+'\\fragment_matches'
            if not os.path.exists(peptide_rep_output_folder):
                os.makedirs(peptide_rep_output_folder)
            fragment_match_dfs = []
            fragment_temp_cutoff = fragment_error_cutoff*1000000
            
            for b in range(0,len(ion_report)):
    
                ion_report_filtered = ion_report.iloc[[b]]
    
                frag_merge_match2 = pd.merge_asof(ion_report_filtered,fragments_extract, left_on='Fragment theoretical monoisotopic mass', 
                                            right_on='Fragment actual monoisotopic mass',
                                            tolerance = fragment_temp_cutoff, allow_exact_matches=True,direction='forward') 
    
                frag_merge_match3 = pd.merge_asof(ion_report_filtered,fragments_extract, left_on='Fragment theoretical monoisotopic mass', 
                                            right_on='Fragment actual monoisotopic mass',
                                            tolerance = fragment_temp_cutoff, allow_exact_matches=True,direction='backward') 
    
                fragment_match_dfs.append(frag_merge_match2)
                fragment_match_dfs.append(frag_merge_match3)  
            
            frag_merge_match = pd.concat(fragment_match_dfs,ignore_index=True)
            output_path_rep = peptide_report_output + '\\frag_merge_match2.csv'
    
            frag_merge_match = frag_merge_match[frag_merge_match['Sequence'].notna()]
            
            # output_path_rep = peptide_report_output + '\\frag_merge_match1.csv'
            # with open(output_path_rep,'w',newline='') as filec:
            #         writerc = csv.writer(filec)
            #         frag_merge_match.to_csv(filec,index=False)
            
            frag_merge_match['Fragment error (Da)'] = abs(frag_merge_match['Fragment actual monoisotopic mass'] - frag_merge_match['Fragment theoretical monoisotopic mass'])
            
            # output_path_rep = peptide_report_output + '\\frag_merge_match_filtered_noerr.csv'
            # with open(output_path_rep,'w',newline='') as filec:
            #         writerc = csv.writer(filec)
            #         frag_merge_match.to_csv(filec,index=False)
                    
            frag_merge_match_filtered = frag_merge_match
    
            
            # output_path_rep = peptide_report_output + '\\frag_merge_match_filtered.csv'
            # with open(output_path_rep,'w',newline='') as filec:
            #         writerc = csv.writer(filec)
            #         frag_merge_match_filtered.to_csv(filec,index=False)
                    
            list_of_df = [g for _, g in frag_merge_match_filtered.groupby(['Scan'])]
            for df in list_of_df:
    
                scan_to_report = int(df['Scan'].values[0])
                peptide = df['Sequence'].values[0]
                seq_cov = seq_coverage_calc(df,ion_report,scan_to_report,peptide)
                unfiltered_psms.append(seq_cov)
    
                issue_seq1 = 'Q(Gln->pyro-Glu)'
                issue_seq2 = 'E(Glu->pyro-Glu)'
                peptide_formal = peptide
                
                if issue_seq1 in peptide_formal:
                    peptide_formal = peptide_formal.replace(issue_seq1,'Q(pyroGlu)')
                elif issue_seq2 in peptide_formal:
                    peptide_formal = peptide_formal.replace(issue_seq2,'E(pyroGlu)')
    
                output_path_rep = peptide_rep_output_folder + '\\' + peptide_formal + '_' + str(scan_to_report) + '_fragment_report.csv'
        
                with open(output_path_rep,'w',newline='') as filec:
                        writerc = csv.writer(filec)
                        df.to_csv(filec,index=False)
        
        #End fragment AMM
        
        e = MSExperiment()
        
        print('here1')
        MzMLFile().load(mzml_path_input, e)
        print('here2')
        unfiltered_psms_df = pd.concat(unfiltered_psms,ignore_index=True)
        print('here3')
        unfiltered_psm_scan_only = unfiltered_psms_df.drop_duplicates(subset='Scan')
        print('here4')
        correlation_results = []
        print('here5')
        for scan in range(0,len(unfiltered_psm_scan_only)):
    
            scan_isolate = unfiltered_psm_scan_only.iloc[[scan]]
            scan_value = scan_isolate['Scan'].values[0]
            scan_filtered_full_results = unfiltered_psms_df[unfiltered_psms_df['Scan'] == scan_value]
    
            for peptide in range(0,len(scan_filtered_full_results)):
                peptide_isolate = scan_filtered_full_results.iloc[[peptide]]
                peptide_value = peptide_isolate['Sequence'].values[0]
                correlation_value = hyperscore_calc(scan_value,peptide_value)
                peptide_isolate['Correlation value'] = correlation_value
                correlation_results.append(peptide_isolate)
                
        all_correlation_results = pd.concat(correlation_results,ignore_index=True)   
        
        output_path_rep = peptide_report_output + '\\all_correlation_results.csv'
        with open(output_path_rep,'w',newline='') as filec:
                writerc = csv.writer(filec)
                all_correlation_results.to_csv(filec,index=False)
        
        return all_correlation_results
        
        subject = 'Your code for '+sample_name+' has finished running'
        text = 'Your code for '+sample_name+' has finished running'
        content = 'Subject: %s\n\n%s' % (subject, text)
        mail = smtplib.SMTP('smtp.gmail.com',587)
        mail.ehlo()
        mail.starttls()
        mail.login('lingjun.li.notifications@gmail.com','eabtnjwaikdssdtd')
        mail.sendmail('lingjun.li.notifications@gmail.com','lawashburn@wisc.edu',content) 
        mail.close()
        
      
    
    subject = 'Your code has finished running'
    text = 'Your code has finished running'
    content = 'Subject: %s\n\n%s' % (subject, text)
    mail = smtplib.SMTP('smtp.gmail.com',587)
    mail.ehlo()
    mail.starttls()
    mail.login('lingjun.li.notifications@gmail.com','eabtnjwaikdssdtd')
    mail.sendmail('lingjun.li.notifications@gmail.com','lawashburn@wisc.edu',content) 
    mail.close()