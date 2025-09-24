# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 11:02:53 2024

@author: lafields2
"""

from pyteomics import mass,parser
import re
import pandas as pd

fasta_to_df = ['LYKPYGW']


variable_mod_dict = {'(DiLeu4plex115)':['K'],
                     '(Sulfo)':['Y']}

# custom_mods = {

max_mods_number = 3

mass_of_H = 1.00794

proton_mass = 1.00727646688
charge = 1 #fragment charge
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
    'K(DiLeu4plex115)' : C*13  + H*27 + N*2 + O*1 + C13*1 + N15*1 + O18*1
    }
termini = {'Standard' : H * 2 + O,
'(Amidated)' : N + H * 3}
PTMs = {'C(Pyro-glu)' : H * -3 - N,
        'Q(Pyro-glu)' : H * -3 - N,
        'E(Pyro-glu)' : H * -3 - N,
        'M(Oxidation)' : O,
        'Y(Sulfo)' : S + O * 3
        }
adducts = {
    'H2O' : H * 2 + O,
    'NH3' : N + H * 3}

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

def check_termini_Key(dict, key):
    if key in dict.keys():
        return dict[key]
    else:
        return termini['Standard']

def check_PTM_Key(dict, key):
    if key in dict.keys():
        return dict[key]

    else:
        print('Key is not in dict')
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
    print('\nnext peptide\n')
    print(pot_mod_peptide)
    number_of_mods = pot_mod_peptide.count('(')
    print(number_of_mods)
    if number_of_mods > 0:
        current_peptide = []
        mass_change_collection = []
        current_peptide.append(pot_mod_peptide)
        print(pot_mod_peptide)
        for a in range(0,number_of_mods):
            peptide_less_mod = current_peptide[-1]
            print('peptide less mod: ',peptide_less_mod)
            ptm_start = (peptide_less_mod.index('('))-1
            ptm_end = (peptide_less_mod.index(')'))+1
            ptm_ID = peptide_less_mod[ptm_start:ptm_end]
            print('PTM ID: ',ptm_ID)
            ptm_mass_change = check_PTM_Key(PTMs, ptm_ID)
            print('PTM mass change: ',ptm_mass_change)
            mass_change_collection.append(ptm_mass_change)
            peptide_less_mod2 = peptide_less_mod[:ptm_start] + peptide_less_mod[ptm_end:]
            print('PTM less mod 2: ',peptide_less_mod2)
            current_peptide.append(peptide_less_mod2)
            
        ptm_mass_change = sum(mass_change_collection)
        
        return ptm_mass_change
    else:
        ptm_mass_change = 0
        return ptm_mass_change

def monoisotopic_mass_calculator(peptide_from_fasta):
        plain_peptide = re.sub("[\(\[].*?[\)\]]", "",peptide_from_fasta)
        res_list_for_fragment = list_of_residues(peptide_from_fasta)
        mass_of_residues = []
        for residue in plain_peptide:
            residue_mass = aa_masses[residue]
            mass_of_residues.append(residue_mass)
        peptide_mass = ((sum(mass_of_residues)) + check_termini(peptide_from_fasta) + check_PTM(peptide_from_fasta))-proton_mass
        mass_to_charge = (peptide_mass + (proton_mass * charge))/charge
        return mass_to_charge

modified_target_sequence = []
for sequence in fasta_to_df:
     forms = parser.isoforms(sequence, variable_mods=variable_mod_dict,max_mods=max_mods_number)   
     forms_list = [*forms]
     for a in forms_list:
         modified_target_sequence.append(a)
         #######
sequences = modified_target_sequence
      #########   

replace_dictionary = {"-Amidated":"(Amidated)", "Ac-":"(Acetylation)","(Glu->pyro-Glu)E":"E(Glu->pyro-Glu)","(Oxidation)M":"M(Oxidation)",
                      "(Gln->pyro-Glu)Q":"Q(Gln->pyro-Glu)","(Sulfo)Y":"Y(Sulfo)","(DiLeu4plex115)K":"K(DiLeu4plex115)"} 

pep_formatted = []
masses = []
for pep in sequences:
    for key, value in replace_dictionary.items():
        pep = pep.replace(key, value)
    pep_formatted.append(pep)

    mass_calc = monoisotopic_mass_calculator(pep)
    masses.append(mass_calc)
    
check = pd.DataFrame()
check['Peptide'] = pep_formatted
check['Mass'] = masses
