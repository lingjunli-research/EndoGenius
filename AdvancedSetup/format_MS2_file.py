# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 14:54:27 2023

@author: lawashburn
"""

import csv
import pandas as pd
import os

# MS2_path = r'C:/Users/lawashburn/Documents/DB_pep_validation/GUI_test_files/2021_0817_CoG_1.ms2'
# output_directory = r"C:\Users\lawashburn\Documents\DB_pep_validation\GUI_test_files\output\v04"

def format_raw_MS2(MS2_path,output_directory):
    backslash_index1 = MS2_path.rfind('\\')
    backslash_index2 = MS2_path.rfind('/')

    if backslash_index1 > backslash_index2:
        backslash_index = backslash_index1
    elif backslash_index2 >= backslash_index1:
        backslash_index = backslash_index2

    base_file_path = MS2_path[0:(backslash_index+1)]

    tissue_type = MS2_path.replace(base_file_path,'')
    tissue_type = tissue_type.replace('.ms2','')
    with open(MS2_path) as input:
        lst = [line.strip() for line in input]
    
    new_list= []
    final_lst = []
    final_lst.append(['m/z', 'resolution', 'charge', 'intensity', 'MS2','scan_number','precursor_charge'])
    ms2_list = []
    
    new = lst
    
    for i in new:
        new_list.append(i.split())
        if '@' in i:
            x = i.split()
            for y in x:
                if '@' in y:
                    ms2 = y[0:y.index('@')]
                    ms2_list.append(str(ms2))
    
    header_list = new_list[0:26]
    new_list = new_list[26:] # starts from line 26 to remove the first few header lines so that program could proceed
    seperation_list = []
    scan_number_list = []    
    precursor_charge_list = []
    
    for i in header_list:
        if 'S' in i:
            scan_number_list.append(i[1])
        if 'Z' in i:
            precursor_charge_list.append(i[1])
    
    for i in range(len(new_list)):
        if 'RetTime' in new_list[i]:
            seperation_list.append(i-1)
        if 'PrecursorInt' in new_list[i]:
            seperation_list.append(i+2)
        if 'S' in new_list[i]:
            scan_number_list.append(new_list[i][1])
        if 'Z' in new_list[i]:
            precursor_charge_list.append(new_list[i][1])
    
    seperation_pairs = []
    start = 0
    for i in range(int(len(seperation_list)/2)):
        seperation_pairs.append((seperation_list[i+start],seperation_list[i+start+1]))
        start +=1 
     
    update_index = 0
    for start,end in seperation_pairs:
        start += update_index
        end += update_index
        new_list[start:end] = '-'
        update_index -= (end-start-1)
    
    ms2_list_index = 0
    scan_number_index = 0
    precursor_charge_index = 0
    
    for element in new_list:
        if element == '-':
            ms2_list_index+=1
            scan_number_index+=1
            precursor_charge_index+=1
            continue   
        element.append(ms2_list[ms2_list_index])
        element.append(scan_number_list[scan_number_index])
        element.append(precursor_charge_list[precursor_charge_index])
        final_lst.append(element)

    out_name = output_directory + '\\'+tissue_type+'_formatted.txt'
    with open(out_name,'w') as output:
        for i in final_lst:
            for j in i:
                output.write(str(j + ','))
            output.write('\n')
    return out_name
# format_raw_MS2(MS2_path,output_directory)
