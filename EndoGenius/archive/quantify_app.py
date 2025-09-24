# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 14:09:44 2024

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



window=Tk()
window.title('Quantification App')
window.geometry('300x300')


single_file = StringVar()

e=Entry(window,show=None,textvariable = single_file,width =50)
e.pack()

def search():
        path_MS2 = askopenfilename(filetypes=[(".MS2 Files",("*.MS2"))]) 
        single_file.set(path_MS2)    

def add():
        var=e.get()
        lb.insert('end',var)
        e.delete(0,'end')

def delete():
        lb.delete(lb.curselection())

def selected_item():
    num_entries = lb.size()
    print(num_entries)
    print(lb.get(0,last=num_entries))

b3=Button(window,text='search',command=search)
b3.pack()

lb=Listbox(window,width=50)
lb.pack()

b1=Button(window,text='add',width=15,command=add)
b1.pack()


b2=Button(window,text='delete',command=delete)
b2.pack()


window.mainloop()
