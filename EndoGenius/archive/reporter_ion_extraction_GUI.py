# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 19:10:47 2024

@author: lafields2
"""

import tkinter as tk
from tkinter import filedialog
from tkinter import StringVar
import pandas as pd
import csv


class QuantificationApp(tk.Tk):
        def __init__(self):
            super().__init__()
    
            self.title("Reporter Ion Extraction")
            self.geometry("760x540")
            self.configure(bg="#423C56")
    
            self.create_widgets()
            
            
    
        def create_widgets(self):
            # Title label
            
            self.exportdirectorypath = StringVar()
            self.importdirectorypath = StringVar()
            self.importspectrapath = StringVar()
            self.errorthreshold = StringVar()
            
            # self.exportdirectorypath.set(r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.1.4")
            # self.importdirectorypath.set(r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.1.2\v10\20240523_DiLeu_TR1_240523184202_5mod_5FDR\20240523_DiLeu_TR1_240523184202\final_results__target.csv")
            # self.importspectrapath.set(r"C:\Users\lawashburn\Documents\EndoGeniusDistributions\version_assessment_output\EndoGenius_v1.1.2\v10\20240523_DiLeu_TR1_240523184202_5mod_5FDR\20240523_DiLeu_TR1_240523184202_formatted.txt")
            # self.errorthreshold.set(20)
            
            title_label = tk.Label(self, text="Reporter Ion Extraction", font=("Inter", 46), bg="#423C56", justify='left',fg="white")
            title_label.pack(pady=0)
    
            frame = tk.Frame(self, bg="#d3d3d3")
            frame.pack(expand=True, fill=tk.BOTH, padx=20, pady=20)
            
            # Import Results
            import_label = tk.Label(frame, text="Import Results:", bg="#d3d3d3", font=("Arial", 12))
            import_label.grid(row=0, column=0, padx=5, pady=5, sticky=tk.E)
    
            self.import_entry = tk.Entry(frame, width=50, font=("Arial", 12),textvariable=self.importdirectorypath)
            self.import_entry.grid(row=0, column=1, padx=5, pady=5)
    
            import_button = tk.Button(frame, text="Browse", command=self.browse_import, font=("Arial", 12))
            import_button.grid(row=0, column=2, padx=5, pady=5)
    
            # Import Results
            import_label = tk.Label(frame, text="Import Spectra:", bg="#d3d3d3", font=("Arial", 12))
            import_label.grid(row=1, column=0, padx=5, pady=5, sticky=tk.E)
    
            self.import_spectra_entry = tk.Entry(frame, width=50, font=("Arial", 12),textvariable=self.importspectrapath)
            self.import_spectra_entry.grid(row=1, column=1, padx=5, pady=5)
    
            import_button = tk.Button(frame, text="Browse", command=self.browse_spectra, font=("Arial", 12))
            import_button.grid(row=1, column=2, padx=5, pady=5)
            
            # IntVar to hold the value of the selected radio button
            self.radio_value = tk.IntVar()
            self.radio_value.set(0)  # Set default value to 0 (no selection)
    
            # Create the radio buttons and associate them with radio_value
            self.radio1 = tk.Radiobutton(frame, text="4-Plex", variable=self.radio_value, value=4,bg="#d3d3d3", font=("Arial", 12),command=self.update_plex,)
            self.radio2 = tk.Radiobutton(frame, text="12-Plex", variable=self.radio_value, value=12,bg="#d3d3d3", font=("Arial", 12),command=self.update_plex,)
            self.radio1.grid(row=2, column=0, padx=5, pady=5)
            self.radio2.grid(row=2, column=1, padx=5, pady=5)
    
            error_label = tk.Label(frame, text="Error tolerance (ppm):", bg="#d3d3d3", font=("Arial", 12))
            error_label.grid(row=2, column=1, padx=0, pady=5, sticky=tk.E)
    
            self.error_entry = tk.Entry(frame, width=10, font=("Arial", 12),textvariable=self.errorthreshold)
            self.error_entry.grid(row=2, column=2, padx=0, pady=5)
    
            self.plex_frame = tk.Frame(frame, bg="#d3d3d3")
            self.plex_frame.grid(row=3, column=0, columnspan=4, padx=5, pady=10)
    
            # Export Directory
            export_label = tk.Label(frame, text="Export Directory:", bg="#d3d3d3", font=("Arial", 12))
            export_label.grid(row=4, column=0, padx=5, pady=5, sticky=tk.E)
    
            self.export_entry = tk.Entry(frame, width=50, font=("Arial", 12),textvariable=self.exportdirectorypath)
            self.export_entry.grid(row=4, column=1, padx=5, pady=5)
    
            export_button = tk.Button(frame, text="Browse", command=self.browse_export, font=("Arial", 12))
            export_button.grid(row=4, column=2, padx=5, pady=5)
    
            analyze_button = tk.Button(frame, text="Begin Analysis", command=self.begin_analysis, font=("Arial", 12))
            analyze_button.grid(row=5, column=0, columnspan=4, pady=10,sticky=tk.N)
    
        def browse_import(self):
            file_path = filedialog.askopenfilename()
            self.import_entry.insert(0, file_path)
        
        def browse_spectra(self):
            file_path = filedialog.askopenfilename()
            self.import_spectra_entry.insert(0, file_path)
    
        def browse_export(self):
            directory_path = filedialog.askdirectory()
            self.export_entry.insert(0, directory_path)
    
        def update_plex(self):
            for widget in self.plex_frame.winfo_children():
                widget.destroy()
    
            plex_count = self.radio_value.get()
            print(plex_count)
            if plex_count == 4:
                
                self.fourplex_1_name_entry = StringVar()
                self.fourplex_1_name_entry.set('115')
                
                self.fourplex_1_mass_entry = StringVar()
                self.fourplex_1_mass_entry.set('115.1253')
                
                self.fourplex_2_name_entry = StringVar()
                self.fourplex_2_name_entry.set('116')
                
                self.fourplex_2_mass_entry = StringVar()
                self.fourplex_2_mass_entry.set('116.1408')
                
                self.fourplex_3_name_entry = StringVar()
                self.fourplex_3_name_entry.set('117')
                
                self.fourplex_3_mass_entry = StringVar()
                self.fourplex_3_mass_entry.set('117.1379')
                
                self.fourplex_4_name_entry = StringVar()
                self.fourplex_4_name_entry.set('118')
                
                self.fourplex_4_mass_entry = StringVar()
                self.fourplex_4_mass_entry.set('118.1534')
                
                plex_name_title_left = tk.Label(self.plex_frame, text="Name", bg="#d3d3d3", font=("Arial", 12))
                plex_name_title_left.grid(row=0, column=1)
                
                plex_mass_title_left = tk.Label(self.plex_frame, text="m/z", bg="#d3d3d3", font=("Arial", 12))
                plex_mass_title_left.grid(row=0, column=2)
                
                plex1_name_label = tk.Label(self.plex_frame, text="1.", bg="#d3d3d3", font=("Arial", 12))
                plex1_name_label.grid(row=1, column=0, padx=5, pady=2, sticky=tk.E)
                
                plex1_name_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.fourplex_1_name_entry)
                plex1_name_entry.grid(row=1, column=1, padx=5, pady=2)
                
                plex1_mass_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.fourplex_1_mass_entry)
                plex1_mass_entry.grid(row=1, column=2, padx=5, pady=2)
                
                plex2_name_label = tk.Label(self.plex_frame, text="2.", bg="#d3d3d3", font=("Arial", 12))
                plex2_name_label.grid(row=2, column=0, padx=5, pady=2, sticky=tk.E)
                
                plex2_name_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.fourplex_2_name_entry)
                plex2_name_entry.grid(row=2, column=1, padx=5, pady=2)
                
                plex2_mass_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.fourplex_2_mass_entry)
                plex2_mass_entry.grid(row=2, column=2, padx=5, pady=2)
                
                plex_name_title_right = tk.Label(self.plex_frame, text="Name", bg="#d3d3d3", font=("Arial", 12))
                plex_name_title_right.grid(row=0, column=4)
                
                plex_mass_title_right = tk.Label(self.plex_frame, text="m/z", bg="#d3d3d3", font=("Arial", 12))
                plex_mass_title_right.grid(row=0, column=5)
                
                plex3_name_label = tk.Label(self.plex_frame, text="3.", bg="#d3d3d3", font=("Arial", 12))
                plex3_name_label.grid(row=1, column=3, padx=5, pady=2, sticky=tk.E)
                
                plex3_name_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.fourplex_3_name_entry)
                plex3_name_entry.grid(row=1, column=4, padx=5, pady=2)
                
                plex3_mass_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.fourplex_3_mass_entry)
                plex3_mass_entry.grid(row=1, column=5, padx=5, pady=2)
                
                plex4_name_label = tk.Label(self.plex_frame, text="4.", bg="#d3d3d3", font=("Arial", 12))
                plex4_name_label.grid(row=2, column=3, padx=5, pady=2, sticky=tk.E)
                
                plex4_name_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.fourplex_4_name_entry)
                plex4_name_entry.grid(row=2, column=4, padx=5, pady=2)
                
                plex4_mass_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.fourplex_4_mass_entry)
                plex4_mass_entry.grid(row=2, column=5, padx=5, pady=2)
                
    
            if plex_count == 12:
                self.twelveplex_1_name_entry = StringVar()
                self.twelveplex_1_name_entry.set('115a')
                
                self.twelveplex_1_mass_entry = StringVar()
                self.twelveplex_1_mass_entry.set('115.12476')
                
                self.twelveplex_2_name_entry = StringVar()
                self.twelveplex_2_name_entry.set('115b')
                
                self.twelveplex_2_mass_entry = StringVar()
                self.twelveplex_2_mass_entry.set('115.13108')
                
                self.twelveplex_3_name_entry = StringVar()
                self.twelveplex_3_name_entry.set('116a')
                
                self.twelveplex_3_mass_entry = StringVar()
                self.twelveplex_3_mass_entry.set('116.12812')
                
                self.twelveplex_4_name_entry = StringVar()
                self.twelveplex_4_name_entry.set('116b')
                
                self.twelveplex_4_mass_entry = StringVar()
                self.twelveplex_4_mass_entry.set('116.13444')
                
                self.twelveplex_5_name_entry = StringVar()
                self.twelveplex_5_name_entry.set('116c')
                
                self.twelveplex_5_mass_entry = StringVar()
                self.twelveplex_5_mass_entry.set('116.14028')
                
                self.twelveplex_6_name_entry = StringVar()
                self.twelveplex_6_name_entry.set('117a')
                
                self.twelveplex_6_mass_entry = StringVar()
                self.twelveplex_6_mass_entry.set('117.13147')
                
                
                
                self.twelveplex_7_name_entry = StringVar()
                self.twelveplex_7_name_entry.set('117b')
                
                self.twelveplex_7_mass_entry = StringVar()
                self.twelveplex_7_mass_entry.set('117.13731')
                
                self.twelveplex_8_name_entry = StringVar()
                self.twelveplex_8_name_entry.set('117c')
                
                self.twelveplex_8_mass_entry = StringVar()
                self.twelveplex_8_mass_entry.set('117.14363')
                
                self.twelveplex_9_name_entry = StringVar()
                self.twelveplex_9_name_entry.set('118a')
                
                self.twelveplex_9_mass_entry = StringVar()
                self.twelveplex_9_mass_entry.set('118.13483')
                
                self.twelveplex_10_name_entry = StringVar()
                self.twelveplex_10_name_entry.set('118b')
                
                self.twelveplex_10_mass_entry = StringVar()
                self.twelveplex_10_mass_entry.set('118.14067')
                
                self.twelveplex_11_name_entry = StringVar()
                self.twelveplex_11_name_entry.set('118c')
                
                self.twelveplex_11_mass_entry = StringVar()
                self.twelveplex_11_mass_entry.set('118.14699')
                
                self.twelveplex_12_name_entry = StringVar()
                self.twelveplex_12_name_entry.set('118d')
                
                self.twelveplex_12_mass_entry = StringVar()
                self.twelveplex_12_mass_entry.set('118.15283')
                
                plex_name_title_left = tk.Label(self.plex_frame, text="Name", bg="#d3d3d3", font=("Arial", 12))
                plex_name_title_left.grid(row=0, column=1)
                
                plex_mass_title_left = tk.Label(self.plex_frame, text="m/z", bg="#d3d3d3", font=("Arial", 12))
                plex_mass_title_left.grid(row=0, column=2)
                
                plex1_name_label = tk.Label(self.plex_frame, text="1.", bg="#d3d3d3", font=("Arial", 12))
                plex1_name_label.grid(row=1, column=0, padx=5, pady=2, sticky=tk.E)
                
                plex1_name_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_1_name_entry)
                plex1_name_entry.grid(row=1, column=1, padx=5, pady=2)
                
                plex1_mass_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_1_mass_entry)
                plex1_mass_entry.grid(row=1, column=2, padx=5, pady=2)
                
                plex2_name_label = tk.Label(self.plex_frame, text="2.", bg="#d3d3d3", font=("Arial", 12))
                plex2_name_label.grid(row=2, column=0, padx=5, pady=2, sticky=tk.E)
                
                plex2_name_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_2_name_entry)
                plex2_name_entry.grid(row=2, column=1, padx=5, pady=2)
                
                plex2_mass_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_2_mass_entry)
                plex2_mass_entry.grid(row=2, column=2, padx=5, pady=2)
                
                plex3_name_label = tk.Label(self.plex_frame, text="3.", bg="#d3d3d3", font=("Arial", 12))
                plex3_name_label.grid(row=3, column=0, padx=5, pady=2, sticky=tk.E)
                
                plex3_name_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_3_name_entry)
                plex3_name_entry.grid(row=3, column=1, padx=5, pady=2)
                
                plex3_mass_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_3_mass_entry)
                plex3_mass_entry.grid(row=3, column=2, padx=5, pady=2)
                
                plex4_name_label = tk.Label(self.plex_frame, text="4.", bg="#d3d3d3", font=("Arial", 12))
                plex4_name_label.grid(row=4, column=0, padx=5, pady=2, sticky=tk.E)
                
                plex4_name_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_4_name_entry)
                plex4_name_entry.grid(row=4, column=1, padx=5, pady=2)
                
                plex4_mass_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_4_mass_entry)
                plex4_mass_entry.grid(row=4, column=2, padx=5, pady=2)
                
                plex5_name_label = tk.Label(self.plex_frame, text="5.", bg="#d3d3d3", font=("Arial", 12))
                plex5_name_label.grid(row=5, column=0, padx=5, pady=2, sticky=tk.E)
                
                plex5_name_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_5_name_entry)
                plex5_name_entry.grid(row=5, column=1, padx=5, pady=2)
                
                plex5_mass_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_5_mass_entry)
                plex5_mass_entry.grid(row=5, column=2, padx=5, pady=2)
                
                plex6_name_label = tk.Label(self.plex_frame, text="6.", bg="#d3d3d3", font=("Arial", 12))
                plex6_name_label.grid(row=6, column=0, padx=5, pady=2, sticky=tk.E)
                
                plex6_name_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_6_name_entry)
                plex6_name_entry.grid(row=6, column=1, padx=5, pady=2)
                
                plex6_mass_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_6_mass_entry)
                plex6_mass_entry.grid(row=6, column=2, padx=5, pady=2)
                
                plex_name_title_right = tk.Label(self.plex_frame, text="Name", bg="#d3d3d3", font=("Arial", 12))
                plex_name_title_right.grid(row=0, column=4)
                
                plex_mass_title_right = tk.Label(self.plex_frame, text="m/z", bg="#d3d3d3", font=("Arial", 12))
                plex_mass_title_right.grid(row=0, column=5)
                
                plex7_name_label = tk.Label(self.plex_frame, text="7.", bg="#d3d3d3", font=("Arial", 12))
                plex7_name_label.grid(row=1, column=3, padx=5, pady=2, sticky=tk.E)
                
                plex7_name_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_7_name_entry)
                plex7_name_entry.grid(row=1, column=4, padx=5, pady=2)
                
                plex7_mass_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_7_mass_entry)
                plex7_mass_entry.grid(row=1, column=5, padx=5, pady=2)
                
                plex8_name_label = tk.Label(self.plex_frame, text="8.", bg="#d3d3d3", font=("Arial", 12))
                plex8_name_label.grid(row=2, column=3, padx=5, pady=2, sticky=tk.E)
                
                plex8_name_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_8_name_entry)
                plex8_name_entry.grid(row=2, column=4, padx=5, pady=2)
                
                plex8_mass_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_8_mass_entry)
                plex8_mass_entry.grid(row=2, column=5, padx=5, pady=2)
                
                plex9_name_label = tk.Label(self.plex_frame, text="9.", bg="#d3d3d3", font=("Arial", 12))
                plex9_name_label.grid(row=3, column=3, padx=5, pady=2, sticky=tk.E)
                
                plex9_name_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_9_name_entry)
                plex9_name_entry.grid(row=3, column=4, padx=5, pady=2)
                
                plex9_mass_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_9_mass_entry)
                plex9_mass_entry.grid(row=3, column=5, padx=5, pady=2)
                
                plex10_name_label = tk.Label(self.plex_frame, text="10.", bg="#d3d3d3", font=("Arial", 12))
                plex10_name_label.grid(row=4, column=3, padx=5, pady=2, sticky=tk.E)
                
                plex10_name_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_10_name_entry)
                plex10_name_entry.grid(row=4, column=4, padx=5, pady=2)
                
                plex10_mass_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_10_mass_entry)
                plex10_mass_entry.grid(row=4, column=5, padx=5, pady=2)
                
                plex11_name_label = tk.Label(self.plex_frame, text="11.", bg="#d3d3d3", font=("Arial", 12))
                plex11_name_label.grid(row=5, column=3, padx=5, pady=2, sticky=tk.E)
                
                plex11_name_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_11_name_entry)
                plex11_name_entry.grid(row=5, column=4, padx=5, pady=2)
                
                plex11_mass_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_11_mass_entry)
                plex11_mass_entry.grid(row=5, column=5, padx=5, pady=2)
                
                plex12_name_label = tk.Label(self.plex_frame, text="12.", bg="#d3d3d3", font=("Arial", 12))
                plex12_name_label.grid(row=6, column=3, padx=5, pady=2, sticky=tk.E)
                
                plex12_name_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_12_name_entry)
                plex12_name_entry.grid(row=6, column=4, padx=5, pady=2)
                
                plex12_mass_entry = tk.Entry(self.plex_frame, width=10, font=("Arial", 12),textvariable=self.twelveplex_12_mass_entry)
                plex12_mass_entry.grid(row=6, column=5, padx=5, pady=2)
    
    
        def begin_analysis(self):
    
            plex_selection = self.radio_value.get()
            tag_name_list = []
            tag_mass_list = []
            
            output_path = self.exportdirectorypath.get()
            results_path = self.importdirectorypath.get()
            spectra_path = self.importspectrapath.get()
            fragment_error_threshold = float(self.errorthreshold.get())
            
            mass_h = 1.00784
            
            if plex_selection == 4:
                ##Retrieve plex names
                tag_name_list.append(self.fourplex_1_name_entry.get())
                tag_name_list.append(self.fourplex_2_name_entry.get())
                tag_name_list.append(self.fourplex_3_name_entry.get())
                tag_name_list.append(self.fourplex_4_name_entry.get())
                ##Retrieve plex masses
                tag_mass_list.append(float(self.fourplex_1_mass_entry.get()))
                tag_mass_list.append(float(self.fourplex_2_mass_entry.get()))
                tag_mass_list.append(float(self.fourplex_3_mass_entry.get()))
                tag_mass_list.append(float(self.fourplex_4_mass_entry.get()))
            
            elif plex_selection == 12:
                ##Retrieve plex names
                tag_name_list.append(self.twelveplex_1_name_entry.get())
                tag_name_list.append(self.twelveplex_2_name_entry.get())
                tag_name_list.append(self.twelveplex_3_name_entry.get())
                tag_name_list.append(self.twelveplex_4_name_entry.get())
                tag_name_list.append(self.twelveplex_5_name_entry.get())
                tag_name_list.append(self.twelveplex_6_name_entry.get())
                tag_name_list.append(self.twelveplex_7_name_entry.get())
                tag_name_list.append(self.twelveplex_8_name_entry.get())
                tag_name_list.append(self.twelveplex_9_name_entry.get())
                tag_name_list.append(self.twelveplex_10_name_entry.get())
                tag_name_list.append(self.twelveplex_11_name_entry.get())
                tag_name_list.append(self.twelveplex_12_name_entry.get())
                ##Retrieve plex masses
                tag_mass_list.append(float(self.twelveplex_1_mass_entry.get()))
                tag_mass_list.append(float(self.twelveplex_2_mass_entry.get()))
                tag_mass_list.append(float(self.twelveplex_3_mass_entry.get()))
                tag_mass_list.append(float(self.twelveplex_4_mass_entry.get()))
                tag_mass_list.append(float(self.twelveplex_5_mass_entry.get()))
                tag_mass_list.append(float(self.twelveplex_6_mass_entry.get()))
                tag_mass_list.append(float(self.twelveplex_7_mass_entry.get()))
                tag_mass_list.append(float(self.twelveplex_8_mass_entry.get()))
                tag_mass_list.append(float(self.twelveplex_9_mass_entry.get()))
                tag_mass_list.append(float(self.twelveplex_10_mass_entry.get()))
                tag_mass_list.append(float(self.twelveplex_11_mass_entry.get()))
                tag_mass_list.append(float(self.twelveplex_12_mass_entry.get()))
                
    
            
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
                    
if __name__ == "__main__":
        app = QuantificationApp()
        app.mainloop()
QuantificationApp()
