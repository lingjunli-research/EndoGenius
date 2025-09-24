# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 16:10:13 2024

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
from tkinter.tix import *
def launch_quant():
    OUTPUT_PATH = Path(__file__).parent
    ASSETS_PATH = OUTPUT_PATH / Path("./assets")
    
    def relative_to_assets(path: str) -> Path:
        return ASSETS_PATH / Path(path)
    
    

    quant_window = Tk()

    quant_window.geometry("760x530")
    quant_window.title('Quantification App')


    quant_menubar = Menu(quant_window)
    quant_filemenu = Menu(quant_menubar, tearoff=0)
    quant_filemenu.add_command(label="Exit", command=quant_window.quit)
    quant_menubar.add_cascade(label="File", menu=quant_filemenu)
    quant_window.configure(bg = "#423C56",menu=quant_menubar)
    quant_output_dir = StringVar()
    quant_input_file_path = StringVar()
    quant_input_sample_name = StringVar()

    def quantify_summary_begin():
        
        quant_output_directory = quant_output_dir.get()

        quant_file_paths = list(quant_entry_4.get(0,END))

        quant_file_names = list(quant_entry_5.get(0,END))

        quant_merge_df = pd.DataFrame()

        if len(quant_file_paths) != len(quant_file_names):
            raise Exception("A file name must be selected for each file input")
            
        else:
            for a in range(0,(len(quant_file_paths))):
                quant_single_file_path = quant_file_paths[a]
                quant_single_name = quant_file_names[a]
                
                quant_single_file = pd.read_csv(quant_single_file_path)
                quant_single_file_filtered = pd.DataFrame()
                quant_single_file_filtered['Peptide'] = quant_single_file['Peptide']
                quant_single_file_filtered[(quant_single_name + ' Intensity')] = quant_single_file['precursor_intensity']
                
                quant_single_file_filtered = quant_single_file_filtered.sort_values(by=(quant_single_name + ' Intensity'),ascending=False)
                quant_single_file_filtered = quant_single_file_filtered.drop_duplicates(subset='Peptide')

                if len(quant_merge_df)==0:
                    quant_merge_df = quant_single_file_filtered
                else:
                    quant_merge_df = pd.merge(quant_merge_df, quant_single_file_filtered, on='Peptide', how='outer')

        quant_file_out_path = quant_output_directory + '\\merged_intensities.csv'
        with open(quant_file_out_path,'w',newline='') as filec:
                writerc = csv.writer(filec)
                quant_merge_df.to_csv(filec,index=False)
        pymsgbox.alert('Your quantification report is complete','Status Update')

    def quant_output_path_get():
        quant_path_out = filedialog.askdirectory() 
        quant_output_dir.set(quant_path_out)

    def quant_input_file_path_get():
        quant_path_input_file = askopenfilename(filetypes=[("CSV Files",("*.csv"))]) 
        quant_input_file_path.set(quant_path_input_file)

    def quant_add_file():
            quant_var=quant_input_file_path.get()
            quant_entry_4.insert('end',quant_var)
            quant_entry_1.delete(0,'end')

    def quant_delete_file():
            quant_entry_4.delete(quant_entry_4.curselection())

    def quant_add_name():
            quant_var=quant_input_sample_name.get()
            quant_entry_5.insert('end',quant_var)
            quant_entry_3.delete(0,'end')

    def quant_delete_sample_name():
            quant_entry_5.delete(quant_entry_5.curselection())


    quant_canvas = Canvas(
        quant_window,
        bg = "#423C56",
        height = 883,
        width = 778,
        bd = 0,
        highlightthickness = 0,
        relief = "ridge"
    )

    quant_canvas.place(x = 0, y = 0)
    quant_canvas.create_rectangle(
        10.0,
        109.0,
        747.0,
        525.0,
        fill="#D9D9D9",
        outline="")

    quant_canvas.create_text(
        15.0,
        128.0,
        anchor="nw",
        text="Import Results:",
        fill="#000000",
        font=("Inter", 16 * -1)
    )

    quant_canvas.create_text(
        5.0,
        0.0,
        anchor="nw",
        text="Quantify Results",
        fill="#FFFFFF",
        font=("Inter", 64 * -1)
    )

    quant_entry_image_1 = PhotoImage(
        file=relative_to_assets("quant_entry_1.png"))
    quant_entry_bg_1 = quant_canvas.create_image(
        301.0,
        138.0,
        image=quant_entry_image_1
    )
    quant_entry_1 = Entry(
        quant_canvas,
        bd=0,
        bg="#FFFFFF",
        highlightthickness=0,
        textvariable=quant_input_file_path
    )
    quant_entry_1.place(
        x=131.0,
        y=123.0,
        width=340.0,
        height=28.0
    )

    quant_button_image_1 = PhotoImage(
        file=relative_to_assets("quant_button_1.png"))
    quant_button_1 = Button(quant_canvas,
        image=quant_button_image_1,
        borderwidth=0,
        highlightthickness=0,
        command=quant_input_file_path_get,
        relief="flat"
    )
    quant_button_1.place(
        x=480.0,
        y=123.0,
        width=78.21710205078125,
        height=30.0
    )

    quant_canvas.create_text(
        11.0,
        452.0,
        anchor="nw",
        text="Export Directory:",
        fill="#000000",
        font=("Inter", 16 * -1)
    )

    quant_entry_image_2 = PhotoImage(
        file=relative_to_assets("quant_entry_2.png"))
    quant_entry_bg_2 = quant_canvas.create_image(
        323.0,
        462.0,
        image=quant_entry_image_2
    )
    quant_entry_2 = Entry(
        quant_canvas,
        bd=0,
        bg="#FFFFFF",
        highlightthickness=0,
        textvariable = quant_output_dir
    )
    quant_entry_2.place(
        x=153.0,
        y=447.0,
        width=340.0,
        height=28.0
    )

    quant_button_image_2 = PhotoImage(
        file=relative_to_assets("quant_button_2.png"))
    quant_button_2 = Button(quant_canvas,
        image=quant_button_image_2,
        borderwidth=0,
        highlightthickness=0,
        command=quant_output_path_get,
        relief="flat"
    )
    quant_button_2.place(
        x=502.0,
        y=447.0,
        width=78.21710205078125,
        height=30.0
    )

    quant_canvas.create_text(
        15.0,
        168.0,
        anchor="nw",
        text="Sample Name:",
        fill="#000000",
        font=("Inter", 16 * -1)
    )

    quant_canvas.create_text(
        14.0,
        200.0,
        anchor="nw",
        text="File List:",
        fill="#000000",
        font=("Inter", 16 * -1)
    )

    quant_entry_image_3 = PhotoImage(
        file=relative_to_assets("quant_entry_3.png"))
    quant_entry_bg_3 = quant_canvas.create_image(
        301.0,
        178.0,
        image=quant_entry_image_3
    )
    quant_entry_3 = Entry(
        quant_canvas,
        bd=0,
        bg="#FFFFFF",
        highlightthickness=0,
        textvariable = quant_input_sample_name
    )
    quant_entry_3.place(
        x=131.0,
        y=163.0,
        width=340.0,
        height=28.0
    )

    quant_entry_4 = Listbox(quant_canvas,width=50)
    quant_entry_4.place(
        x=17.0,
        y=225.0,
        width=340.0,
        height=144.0
    )

    quant_canvas.create_text(
        352.0,
        200.0,
        anchor="nw",
        text="Name List:",
        fill="#000000",
        font=("Inter", 16 * -1)
    )

    quant_entry_5 = Listbox(quant_canvas,width=50)
    quant_entry_5.place(
        x=372.0,
        y=225.0,
        width=340.0,
        height=144.0
    )

    quant_button_image_3 = PhotoImage(
        file=relative_to_assets("quant_button_3.png"))
    quant_button_3 = Button(quant_canvas,
        image=quant_button_image_3,
        borderwidth=0,
        highlightthickness=0,
        command=quant_add_name,
        relief="flat"
    )
    quant_button_3.place(
        x=480.0,
        y=162.0,
        width=109.44155883789062,
        height=31.0
    )

    quant_button_image_4 = PhotoImage(
        file=relative_to_assets("quant_button_4.png"))
    quant_button_4 = Button(quant_canvas,
        image=quant_button_image_4,
        borderwidth=0,
        highlightthickness=0,
        command=quant_add_file,
        relief="flat"
    )
    quant_button_4.place(
        x=558.0,
        y=123.0,
        width=77.21710205078125,
        height=30.0
    )

    quant_button_image_5 = PhotoImage(
        file=relative_to_assets("quant_button_5.png"))
    quant_button_5 = Button(quant_canvas,
        image=quant_button_image_5,
        borderwidth=0,
        highlightthickness=0,
        command=quant_delete_file,
        relief="flat"
    )
    quant_button_5.place(
        x=131.0,
        y=381.0,
        width=100.0,
        height=30.0
    )

    quant_button_image_6 = PhotoImage(
        file=relative_to_assets("quant_button_6.png"))
    quant_button_6 = Button(quant_canvas,
        image=quant_button_image_6,
        borderwidth=0,
        highlightthickness=0,
        command=quant_delete_sample_name,
        relief="flat"
    )
    quant_button_6.place(
        x=503.0,
        y=381.0,
        width=100.0,
        height=30.0
    )

    quant_button_image_7 = PhotoImage(
        file=relative_to_assets("quant_button_7.png"))
    quant_button_7 = Button(quant_canvas,
        image=quant_button_image_7,
        borderwidth=0,
        highlightthickness=0,
        command=quantify_summary_begin,
        relief="flat"
    )
    quant_button_7.place(
        x=319.0,
        y=487.0,
        width=141.0,
        height=30.0
    )

    quant_canvas.create_rectangle(
        638.0,
        828.0,
        738.0,
        928.0,
        fill="#000000",
        outline="")
    quant_window.resizable(False, False)
    quant_window.mainloop()
