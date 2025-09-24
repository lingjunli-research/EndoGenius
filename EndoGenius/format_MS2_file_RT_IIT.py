# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 16:55:58 2024

@author: lawashburn
"""

import csv
import pandas as pd
import os
import pathlib
MS2_path = r"D:\Manuscripts\2024_12plexNP_DiLeu\20240718_480_DiLeu\Mitch_Polyamine_TR3.ms2"
output_directory = r"D:\Manuscripts\2024_12plexNP_DiLeu\20240718_480_DiLeu"
def format_raw_MS2(MS2_path,output_directory):
    
    def remove_lines_starting_with_Z(input_file, output_file):
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                if not line.startswith('Z'):
                    outfile.write(line)
                
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
    final_lst.append(['fragment_mz', 
                      'fragment_intensity', 
                      'fragment_z', 
                      'fragment_resolution', 
                      'precursor_mz',
                      'ms2_scan',
                      'precursor_z',
                      'precursor_RT',
                      'IonInjectTime',
                      'ms1_scan',
                      'precursor_intensity'])
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
    rt_list = []
    iit_list = []
    precursor_scan_list = []
    precursor_int_list = []
    
    for i in header_list:
        
    
        if 'S' in i:
            scan_number_list.append(i[1])
        if 'Z' in i:
            precursor_charge_list.append(i[1])
        if 'RetTime' in i:
            rt_list.append(i[2])
        if 'IonInjectionTime' in i:
            iit_list.append(i[2])
        if 'PrecursorScan' in i:
            precursor_scan_list.append(i[2])
        if 'PrecursorInt' in i:
            precursor_int_list.append(i[2])
    
    for i in range(len(new_list)):
        if 'RetTime' in new_list[i]:
            seperation_list.append(i-1)
        if 'PrecursorInt' in new_list[i]:
            seperation_list.append(i+2)
        if 'S' in new_list[i]:
            scan_number_list.append(new_list[i][1])
        if 'Z' in new_list[i]:
            precursor_charge_list.append(new_list[i][1])
        if 'RetTime' in new_list[i]:
            rt_list.append(new_list[i][2])
        if 'IonInjectionTime' in new_list[i]:
            iit_list.append(new_list[i][2])
        if 'PrecursorScan' in new_list[i]:
            precursor_scan_list.append(new_list[i][2])
        if 'PrecursorInt' in new_list[i]:
            precursor_int_list.append(new_list[i][2])
    
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
    rt_index = 0
    iit_index = 0
    precursor_scan_index = 0
    precursor_intensity_index = 0
    
    for element in new_list:
        if element == '-':
            ms2_list_index+=1
            scan_number_index+=1
            precursor_charge_index+=1
            rt_index+=1
            iit_index+=1
            precursor_scan_index+=1
            precursor_intensity_index+=1
            continue   
        element.append(ms2_list[ms2_list_index])
        element.append(scan_number_list[scan_number_index])
        element.append(precursor_charge_list[precursor_charge_index])
        element.append(rt_list[rt_index])
        element.append(iit_list[iit_index])
        element.append(precursor_scan_list[precursor_scan_index])
        element.append(precursor_int_list[precursor_intensity_index])
        final_lst.append(element)
    out_name = output_directory + '\\'+tissue_type+'_formatted_draft.txt'

    with open(out_name,'w') as output:
        for i in final_lst:
            for j in i:
                output.write(str(j + ','))
            output.write('\n')
    
    input_file = out_name  # Replace with your input file path
    output_file = output_directory + '\\'+tissue_type+'_formatted.txt'  # Replace with your desired output file path
    remove_lines_starting_with_Z(input_file, output_file)
    
    return output_file

format_raw_MS2(MS2_path,output_directory)