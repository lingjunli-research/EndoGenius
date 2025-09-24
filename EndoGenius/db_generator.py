# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 17:01:10 2023

@author: lawashburn
"""

from Bio.SeqIO.FastaIO import SimpleFastaParser
import random
from pyteomics import mass
import re
import pandas as pd
import csv
from pyteomics import parser

def make_a_DB(variable_mod_dict,fasta_path,output_folder,max_mods_number):


    decoy_algorithm = 'Shuffle'

    mass_of_H = 1.00794
    
    proton_mass = 1.00727646688
    charge = 1  # fragment charge
    H = 1.0078250352
    O = 15.99491463
    C = 12.0000000
    N = 14.003074
    P = 30.973762
    S = 31.9720707
    C13 = 13.003355
    O18 = 17.9991610 
    N15 = 15.0001
    H2 = 2.014101777844
    X = 0
    
    replace_dictionary = {
        "-Amidated": "(Amidated)", 
        "Ac-": "(Acetylation)",
        "(Pyro-glu)E": "E(Pyro-glu)",
        "(Gln->pyro-Glu)Q": "Q(Pyro-glu)",
        "(Glu->pyro-Glu)E": "E(Pyro-glu)",
        "(Oxidation)M": "M(Oxidation)",
        "(Pyro-glu)Q": "Q(Pyro-glu)",
        "(Sulfo)Y": "Y(Sulfo)",
        "(Phospho)S": "S(Phospho)",
        "(12PlexDiLeu)K": "K(12PlexDiLeu)",  # Ensure the custom modification is replaced correctly
        "12PlexDiLeu-": "(12PlexDiLeu)"  # Ensure the custom N-terminal modification is replaced correctly
    }
    
    
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
        'Q(Pyro-glu)' : C*5  + H*5  + N + O*2,
        'E(Pyro-glu)' : C*5  + H*4 + O*3,
        'M(Oxidation)' : C*5  + H*9  + N   + O*2   + S,
        'Y(Sulfo)' :  C*9  + H*9  + N   + O*5 + S,
        'S(Phospho)' : C*3  + H*6  + N   + O*5 + P,
        'X': X,
        'K(12PlexDiLeu)': H*15 + C*7 + C13*1 + N15*1 + O18*1  # Add custom modification here
        }
    n_term_mods = {
        'Standard': X,
        '(12PlexDiLeu)': H*15 + C*7 + C13 + N15 + O18  # Add custom N-terminal modification here
    }

    c_term_mods = {
        'Standard': O + H*2,
        '(Amidated)': N + H * 3
    }
    
    PTMs = {
        'C(Pyro-glu)': H * -3 - N,
        'Q(Pyro-glu)': H * -3 - N,
        'E(Pyro-glu)': H * -3 - N,
        'M(Oxidation)': O,
        'Y(Sulfo)': S + O * 3,
        'S(Phospho)' : C*3  + H*6  + N   + O*5 + P,
        'K(12PlexDiLeu)': H*15 + C*7 + C13*1 + N15*1 + O18*1  # Add custom modification here as well
    }
    adducts = {
        'H2O' : H * 2 + O,
        'NH3' : N + H * 3}
    
    def check_term_Key(term_dict, key):
        if key in term_dict.keys():
            return term_dict[key]
        else:
            return term_dict['Standard']

    def check_n_term(pot_mod_peptide):
        if pot_mod_peptide.startswith('('):
            term_end = pot_mod_peptide.index(')') + 1
            n_term_ID = pot_mod_peptide[:term_end]
            n_term_mass_change = check_term_Key(n_term_mods, n_term_ID)
            return n_term_mass_change
        else:
            return n_term_mods['Standard']

    def check_c_term(pot_mod_peptide):
        if pot_mod_peptide.endswith(')'):
            term_start = pot_mod_peptide.rfind('(')
            c_term_ID = pot_mod_peptide[term_start:]
            c_term_mass_change = check_term_Key(c_term_mods, c_term_ID)
            return c_term_mass_change
        else:
            return c_term_mods['Standard']
    
    def check_PTM_Key(dict, key):
        if key in dict.keys():
            return dict[key]
        else:
            return 0

    def check_PTM(pot_mod_peptide):
        mod_masses = []
        pattern = re.compile(r'([A-Z])\(([^)]+)\)')
        for match in pattern.finditer(pot_mod_peptide):
            res, mod = match.groups()
            ptm_mass_change = check_PTM_Key(PTMs, f"{res}({mod})")
            mod_masses.append(ptm_mass_change)
        return sum(mod_masses)
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
    def monoisotopic_mass_calculator(peptide_from_fasta):
        plain_peptide = re.sub("[\(\[].*?[\)\]]", "", peptide_from_fasta)
        res_list_for_fragment = list_of_residues(peptide_from_fasta)
        mass_of_residues = []
        for residue in plain_peptide:
            residue_mass = aa_masses[residue]
            mass_of_residues.append(residue_mass)
        peptide_mass = ((sum(mass_of_residues)) + check_n_term(peptide_from_fasta) + check_c_term(peptide_from_fasta) + check_PTM(peptide_from_fasta)) - proton_mass
        mass_to_charge = (peptide_mass + (proton_mass * charge)) / charge
        return mass_to_charge
    
    fasta_to_df = []
    with open(fasta_path) as fasta_file:  # Will close handle cleanly
        for title, sequence in SimpleFastaParser(fasta_file):
            fasta_to_df.append(sequence)
     
    modified_target_sequence = []
    for sequence in fasta_to_df:
        forms = parser.isoforms(sequence, variable_mods=variable_mod_dict, max_mods=max_mods_number)   
        forms_list = [*forms]
        for a in forms_list:
            modified_target_sequence.append(a)

    decoy_df_storage = []        
    if decoy_algorithm == 'Shuffle':
        for seq in fasta_to_df:
            seq_shuffle = ''.join(random.sample(seq, len(seq)))
            if seq_shuffle not in decoy_df_storage:
                decoy_df_storage.append(seq_shuffle)
            else:
                seq_shuffle = ''.join(random.sample(seq, len(seq)))
                decoy_df_storage.append(seq_shuffle)
    
    if decoy_algorithm == 'Reverse':
        for seq in fasta_to_df:
            seq_reverse = seq[::-1]
            decoy_df_storage.append(seq_reverse)
    
    modified_decoy_sequence = []
    for sequence in decoy_df_storage:
         forms = parser.isoforms(sequence, variable_mods=variable_mod_dict,max_mods=max_mods_number)   
         forms_list = [*forms]
         for a in forms_list:
             modified_decoy_sequence.append(a)
    
    sequences = modified_target_sequence + modified_decoy_sequence
    
    
    pep_formatted = []
    masses = []
    for pep in sequences:
        for key, value in replace_dictionary.items():
            pep = pep.replace(key, value)
        pep_formatted.append(pep)
        
        mass_calc = monoisotopic_mass_calculator(pep)
        masses.append(mass_calc)
    
    df = pd.DataFrame()
    df['Sequence'] = pep_formatted
    df['Precursor theoretical monoisotopic mass'] = masses
    
    df = df[~df['Sequence'].str.contains(r'\)\(')]
    
    file_path = output_folder + '\\'+ decoy_algorithm + '_database.csv'
    with open(file_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            df.to_csv(filec,index=False) 
            
    return file_path
