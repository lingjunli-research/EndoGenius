# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 14:42:35 2024

@author: lafields2
"""

from pathlib import Path
from tkinter.filedialog import askopenfilename
from tkinter import filedialog
import pymsgbox
# from tkinter import *
# Explicit imports to satisfy Flake8
from tkinter import Tk, Canvas, Entry, Text, Button, PhotoImage, Checkbutton, messagebox
from tkinter import *
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
import csv
import webbrowser
import os
#from tkinter.tix import *

import subprocess

def reporter_ion_extraction_begin():
    extract_window = Tk()
    extract_window.geometry("750x670")
    extract_window.configure(bg = "#423C56")
    extract_window.attributes("-topmost", True)
    extract_window.title('Reporter Ion Extraction')
    
    exportdirectorypath = StringVar()
    importdirectorypath = StringVar()
    importspectrapath = StringVar()
    errorthreshold = StringVar()
    
    exportdirectorypath.set(r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.1.4")
    importdirectorypath.set(r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.1.2\v10\20240523_DiLeu_TR1_240523184202_5mod_5FDR\20240523_DiLeu_TR1_240523184202\final_results__target.csv")
    importspectrapath.set(r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.1.2\v10\20240523_DiLeu_TR1_240523184202_5mod_5FDR\20240523_DiLeu_TR1_240523184202_formatted.txt")
    errorthreshold.set(20)
    
    twelveplex_1_name_entry = StringVar()
    twelveplex_2_name_entry = StringVar()
    twelveplex_3_name_entry = StringVar()
    twelveplex_4_name_entry = StringVar()
    twelveplex_5_name_entry = StringVar()
    twelveplex_6_name_entry = StringVar()
    twelveplex_7_name_entry = StringVar()
    twelveplex_8_name_entry = StringVar()
    twelveplex_9_name_entry = StringVar()
    twelveplex_10_name_entry = StringVar()
    twelveplex_11_name_entry = StringVar()
    twelveplex_12_name_entry = StringVar()
    
    twelveplex_1_mass_entry = StringVar()
    twelveplex_2_mass_entry = StringVar()
    twelveplex_3_mass_entry = StringVar()
    twelveplex_4_mass_entry = StringVar()
    twelveplex_5_mass_entry = StringVar()
    twelveplex_6_mass_entry = StringVar()
    twelveplex_7_mass_entry = StringVar()
    twelveplex_8_mass_entry = StringVar()
    twelveplex_9_mass_entry = StringVar()
    twelveplex_10_mass_entry = StringVar()
    twelveplex_11_mass_entry = StringVar()
    twelveplex_12_mass_entry = StringVar()
    
    twelveplex_1_name_entry.set('115a')
    twelveplex_2_name_entry.set('115b')
    twelveplex_3_name_entry.set('116a')
    twelveplex_4_name_entry.set('116b')
    twelveplex_5_name_entry.set('116c')
    twelveplex_6_name_entry.set('117a')
    twelveplex_7_name_entry.set('117b')
    twelveplex_8_name_entry.set('117c')
    twelveplex_9_name_entry.set('118a')
    twelveplex_10_name_entry.set('118b')
    twelveplex_11_name_entry.set('118c')
    twelveplex_12_name_entry.set('118d')
    
    twelveplex_1_mass_entry.set('115.12476')
    twelveplex_2_mass_entry.set('115.13108')
    twelveplex_3_mass_entry.set('116.12812')
    twelveplex_4_mass_entry.set('116.13444')
    twelveplex_5_mass_entry.set('116.14028')
    twelveplex_6_mass_entry.set('117.13147')
    twelveplex_7_mass_entry.set('117.13731')
    twelveplex_8_mass_entry.set('117.14363')
    twelveplex_9_mass_entry.set('118.13483')
    twelveplex_10_mass_entry.set('118.14067')
    twelveplex_11_mass_entry.set('118.14699')
    twelveplex_12_mass_entry.set('118.15283')
    
    canvas = Canvas(extract_window,bg = "#423C56",height = 800,width = 729,bd = 0,highlightthickness = 0,relief = "ridge")
        
    canvas.place(x = 0, y = 0)
    canvas.create_text(19.0,9.0,anchor="nw",text="Reporter Ion Extraction",fill="#FFFFFF",font=("Inter", 64 * -1))
    
    x=19
    y=90
    width= 700
    height = 550
    
    canvas.create_rectangle(x, y, x+width, y+height,fill="#D9D9D9",outline="")
    
    def eg_results_path():
        eg_results_path_csv = askopenfilename(filetypes=[("CSV Files",("*.csv"))]) 
        importdirectorypath.set(eg_results_path_csv)
    
    def eg_results_spectra_path():
        eg_results_path_txt = askopenfilename(filetypes=[("Text Files",("*.txt"))]) 
        importspectrapath.set(eg_results_path_txt)
        
    def export_folder():
        filename = filedialog.askdirectory()
        exportdirectorypath.set(filename)
    
    def begin_reporter_ion_extraction():
        tag_name_list = []
        tag_mass_list = []
        
        output_path = exportdirectorypath.get()
        results_path = importdirectorypath.get()
        spectra_path = importspectrapath.get()
        fragment_error_threshold = float(errorthreshold.get())
        
        mass_h = 1.00784
        
        tag_name_list.append(twelveplex_1_name_entry.get())
        tag_name_list.append(twelveplex_2_name_entry.get())
        tag_name_list.append(twelveplex_3_name_entry.get())
        tag_name_list.append(twelveplex_4_name_entry.get())
        tag_name_list.append(twelveplex_5_name_entry.get())
        tag_name_list.append(twelveplex_6_name_entry.get())
        tag_name_list.append(twelveplex_7_name_entry.get())
        tag_name_list.append(twelveplex_8_name_entry.get())
        tag_name_list.append(twelveplex_9_name_entry.get())
        tag_name_list.append(twelveplex_10_name_entry.get())
        tag_name_list.append(twelveplex_11_name_entry.get())
        tag_name_list.append(twelveplex_12_name_entry.get())
        ##Retrieve plex masses
        tag_mass_list.append(float(twelveplex_1_mass_entry.get()))
        tag_mass_list.append(float(twelveplex_2_mass_entry.get()))
        tag_mass_list.append(float(twelveplex_3_mass_entry.get()))
        tag_mass_list.append(float(twelveplex_4_mass_entry.get()))
        tag_mass_list.append(float(twelveplex_5_mass_entry.get()))
        tag_mass_list.append(float(twelveplex_6_mass_entry.get()))
        tag_mass_list.append(float(twelveplex_7_mass_entry.get()))
        tag_mass_list.append(float(twelveplex_8_mass_entry.get()))
        tag_mass_list.append(float(twelveplex_9_mass_entry.get()))
        tag_mass_list.append(float(twelveplex_10_mass_entry.get()))
        tag_mass_list.append(float(twelveplex_11_mass_entry.get()))
        tag_mass_list.append(float(twelveplex_12_mass_entry.get()))
        
        spectra = pd.read_csv(spectra_path, sep=",",skiprows=[0], names= ['fragment_mz',
                                                                                      'fragment_intensity',
                                                                                      'fragment_z',
                                                                                      'fragment_resolution',
                                                                                      'precursor_mz',
                                                                                      'ms2_scan',
                                                                                      'precursor_z',
                                                                                      'precursor_RT',
                                                                                      'IonInjectTime',
                                                                                      'ms1_scan',
                                                                                      'precursor_intensity',
                                                                                      'null'])
    
        spectra['fragment_z'] = spectra['fragment_z'].replace(0, 1)
        spectra['fragment_monoisotopic_mass'] = (spectra['fragment_mz'] * spectra['fragment_z']) - (mass_h * spectra['fragment_z'])
    
        results = pd.read_csv(results_path)
        
        spectra_w_details = pd.DataFrame()
    
        peptide_name_log = []
        tag_name_log = []
        tag_intensity_log = []
        scan_log = []
        fragment_mz_log = []
    
        for y in tag_mass_list:
            y_monoisotopic = (y*1) - (mass_h*1)
            spectra[str(y) + ' error'] = ((abs(spectra['fragment_monoisotopic_mass'] - y_monoisotopic))/y_monoisotopic)*1E6
            spectra_w_details = spectra
    
        for a in range(0,len(results)):
            peptide_id = results['Peptide'].iloc[a]
            if '(12PlexDiLeu)' in peptide_id:
                scan_id = results['Scan'].iloc[a]
    
                spectra_filtered = spectra_w_details[spectra_w_details['ms2_scan'] == scan_id]
                
                for k in range(0,len(tag_name_list)):
                    tag_name_selected = tag_name_list[k]
                    tag_mass_selected = tag_mass_list[k]
                    
                    if( spectra_filtered[str(tag_mass_selected) + ' error'].min()) <= fragment_error_threshold:
                        spectra_filtered_tag = spectra_filtered[spectra_filtered[str(tag_mass_selected) + ' error'] <= fragment_error_threshold]
                        tag_intensity_log.append(spectra_filtered_tag['fragment_intensity'].max())
                        tag_name_log.append(tag_name_selected)
                        peptide_name_log.append(peptide_id)
                        scan_log.append(scan_id)
                        
                        filter_filter_df = spectra_filtered_tag[spectra_filtered_tag['fragment_intensity'] == (spectra_filtered_tag['fragment_intensity'].max())]
                        fragment_mz_log.append(filter_filter_df['fragment_mz'].iloc[0])
                        
                    else:
                        tag_intensity_log.append(0)
                        tag_name_log.append(tag_name_selected)
                        peptide_name_log.append(peptide_id)
                        scan_log.append(scan_id)
                        fragment_mz_log.append(0)
                
        report_ion_intensity_df = pd.DataFrame()
        report_ion_intensity_df['Peptide'] = peptide_name_log
        report_ion_intensity_df['Scan'] = scan_log
        report_ion_intensity_df['Tag'] = tag_name_log
        report_ion_intensity_df['Intensity'] = tag_intensity_log
        report_ion_intensity_df['Fragment_mz'] = fragment_mz_log
    
        # Pivot the table without hardcoding the values in the 'Tag' column
        pivot_df = report_ion_intensity_df.pivot_table(index=['Peptide', 'Scan'], columns='Tag', values='Intensity').reset_index()
    
        # Flatten the MultiIndex columns
        pivot_df.columns.name = None
        pivot_df.columns = ['Peptide', 'Scan'] + [f'Intensity_{tag}' for tag in pivot_df.columns[2:]]
    
    
        output_path_rep = output_path + '\\reporter_ions_extracted.csv'
        # df.reset_index().to_feather(output_path_rep) 
    
        with open(output_path_rep,'w',newline='') as filec:
                writerc = csv.writer(filec)
                pivot_df.to_csv(filec,index=False)       
        
        messagebox.showinfo("Process Complete", "Reporter ions have been extracted.")
        
    canvas.create_text(26.0,100.0,anchor="nw",text="EndoGenius Results Path: ",fill="#000000",font=("Inter", 16 * -1))
    eg_results_entry = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=importdirectorypath)
    eg_results_entry.place(x=274.0,y=100.0,width=340.0,height=28.0)
    eg_results_button_1 = Button(canvas,text='Browse',borderwidth=0,highlightthickness=0,relief="flat",command=eg_results_path)
    eg_results_button_1.place(x=625.0,y=100.0,width=77.21710205078125,height=30.0)
    
    canvas.create_text(110.0,150.0,anchor="nw",text="Spectra Path: ",fill="#000000",font=("Inter", 16 * -1))
    eg_spectra_entry = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=importspectrapath)
    eg_spectra_entry.place(x=274.0,y=150.0,width=340.0,height=28.0)
    eg_spectra_button_1 = Button(canvas,text='Browse',borderwidth=0,highlightthickness=0,relief="flat",command=eg_results_spectra_path)
    eg_spectra_button_1.place(x=625.0,y=150.0,width=77.21710205078125,height=30.0)
    
    canvas.create_text(115.0,205.0,anchor="nw",text="Export Directory: ",fill="#000000",font=("Inter", 16 * -1))
    eg_export_entry = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=exportdirectorypath)
    eg_export_entry.place(x=274.0,y=200.0,width=340.0,height=28.0)
    eg_export_button_1 = Button(canvas,text='Browse',borderwidth=0,highlightthickness=0,relief="flat",command=export_folder)
    eg_export_button_1.place(x=625.0,y=200.0,width=77.21710205078125,height=30.0)
    
    canvas.create_text(110.0,255.0,anchor="nw",text="Error tolerance (ppm): ",fill="#000000",font=("Inter", 16 * -1))
    error_entry = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=errorthreshold)
    error_entry.place(x=274.0,y=250.0,width=100.0,height=28.0)
    
    canvas.create_text(75.0,290.0,anchor="nw",text="Name",fill="#000000",font=("Inter", 16 * -1))
    canvas.create_text(215.0,290.0,anchor="nw",text="m/z",fill="#000000",font=("Inter", 16 * -1))
    
    canvas.create_text(30.0,314.0,anchor="nw",text="1: ",fill="#000000",font=("Inter", 16 * -1))
    entry115a_name = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_1_name_entry)
    entry115a_name.place(x=50.0,y=310.0,width=100.0,height=28.0)
    entry115a_mass = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_1_mass_entry)
    entry115a_mass.place(x=175.0,y=310.0,width=100.0,height=28.0)
    
    canvas.create_text(30.0,355.0,anchor="nw",text="2: ",fill="#000000",font=("Inter", 16 * -1))
    entry115b_name = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_2_name_entry)
    entry115b_name.place(x=50.0,y=350.0,width=100.0,height=28.0)
    entry115b_mass = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_2_mass_entry)
    entry115b_mass.place(x=175.0,y=350.0,width=100.0,height=28.0)
    
    canvas.create_text(30.0,393.0,anchor="nw",text="3: ",fill="#000000",font=("Inter", 16 * -1))
    entry116a_name = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_3_name_entry)
    entry116a_name.place(x=50.0,y=390.0,width=100.0,height=28.0)
    entry116a_mass = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_3_mass_entry)
    entry116a_mass.place(x=175.0,y=390.0,width=100.0,height=28.0)
    
    canvas.create_text(30.0,435.0,anchor="nw",text="4: ",fill="#000000",font=("Inter", 16 * -1))
    entry116b_name = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_4_name_entry)
    entry116b_name.place(x=50.0,y=430.0,width=100.0,height=28.0)
    entry116b_mass = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_4_mass_entry)
    entry116b_mass.place(x=175.0,y=430.0,width=100.0,height=28.0)
    
    canvas.create_text(30.0,475.0,anchor="nw",text="5: ",fill="#000000",font=("Inter", 16 * -1))
    entry116c_name = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_5_name_entry)
    entry116c_name.place(x=50.0,y=470.0,width=100.0,height=28.0)
    entry116c_mass = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_5_mass_entry)
    entry116c_mass.place(x=175.0,y=470.0,width=100.0,height=28.0)
    
    canvas.create_text(30.0,515.0,anchor="nw",text="6: ",fill="#000000",font=("Inter", 16 * -1))
    entry116d_name = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_6_name_entry)
    entry116d_name.place(x=50.0,y=510.0,width=100.0,height=28.0)
    entry116d_mass = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_6_mass_entry)
    entry116d_mass.place(x=175.0,y=510.0,width=100.0,height=28.0)
    
    canvas.create_text(400.0,290.0,anchor="nw",text="Name",fill="#000000",font=("Inter", 16 * -1))
    canvas.create_text(540.0,290.0,anchor="nw",text="m/z",fill="#000000",font=("Inter", 16 * -1))
    
    canvas.create_text(355.0,314.0,anchor="nw",text="7: ",fill="#000000",font=("Inter", 16 * -1))
    entry115a_name = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_7_name_entry)
    entry115a_name.place(x=375.0,y=310.0,width=100.0,height=28.0)
    entry115a_mass = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_7_mass_entry)
    entry115a_mass.place(x=500.0,y=310.0,width=100.0,height=28.0)
    
    canvas.create_text(355.0,355.0,anchor="nw",text="8: ",fill="#000000",font=("Inter", 16 * -1))
    entry115b_name = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_8_name_entry)
    entry115b_name.place(x=375.0,y=350.0,width=100.0,height=28.0)
    entry115b_mass = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_8_mass_entry)
    entry115b_mass.place(x=500.0,y=350.0,width=100.0,height=28.0)
    
    canvas.create_text(355.0,393.0,anchor="nw",text="9: ",fill="#000000",font=("Inter", 16 * -1))
    entry116a_name = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_9_name_entry)
    entry116a_name.place(x=375.0,y=390.0,width=100.0,height=28.0)
    entry116a_mass = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_9_mass_entry)
    entry116a_mass.place(x=500.0,y=390.0,width=100.0,height=28.0)
    
    canvas.create_text(347.0,435.0,anchor="nw",text="10: ",fill="#000000",font=("Inter", 16 * -1))
    entry116b_name = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_10_name_entry)
    entry116b_name.place(x=375.0,y=430.0,width=100.0,height=28.0)
    entry116b_mass = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_10_mass_entry)
    entry116b_mass.place(x=500.0,y=430.0,width=100.0,height=28.0)
    
    canvas.create_text(347.0,475.0,anchor="nw",text="11: ",fill="#000000",font=("Inter", 16 * -1))
    entry116c_name = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_11_name_entry)
    entry116c_name.place(x=375.0,y=470.0,width=100.0,height=28.0)
    entry116c_mass = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_11_mass_entry)
    entry116c_mass.place(x=500.0,y=470.0,width=100.0,height=28.0)
    
    canvas.create_text(347.0,515.0,anchor="nw",text="12: ",fill="#000000",font=("Inter", 16 * -1))
    entry116d_name = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_12_name_entry)
    entry116d_name.place(x=375.0,y=510.0,width=100.0,height=28.0)
    entry116d_mass = Entry(canvas,bd=0,bg="#FFFFFF",highlightthickness=0,textvariable=twelveplex_12_mass_entry)
    entry116d_mass.place(x=500.0,y=510.0,width=100.0,height=28.0)
    
    button_2 = Button(canvas,text='Begin Analysis',borderwidth=0,highlightthickness=0,relief="flat",command=begin_reporter_ion_extraction)
    button_2.place(x=300.0,y=575.0,width=90,height=30.0)
    
    extract_window.resizable(True, True)
    extract_window.mainloop()
    
reporter_ion_extraction_begin()
