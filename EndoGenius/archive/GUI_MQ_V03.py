# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 16:49:56 2024

@author: lafields2
"""

from pathlib import Path
# Explicit imports to satisfy Flake8
from tkinter import Tk, Canvas, Entry, Button, PhotoImage, filedialog
from tkinter import *

OUTPUT_PATH = Path(__file__).parent
ASSETS_PATH = OUTPUT_PATH / Path("./assets")


def relative_to_assets(path: str) -> Path:
    return ASSETS_PATH / Path(path)




window = Tk()
window.geometry("800x580")
window.configure(bg = "#423C56")
window.title('MotifQuest')

min_motif_len_in = StringVar()
min_motif_len_in.set('3')

min_num_motif_inst = StringVar()
min_num_motif_inst.set('2')

part_motif_flank = StringVar()
part_motif_flank.set('3')

fasta_path_input = StringVar()
fasta_path_input.set(r"D:\Manuscripts\2024_MotifQuest\MotifQuest_for_EG\input\duplicate_removed_crustacean_database_validated_formatted20220725.fasta")

t_val= StringVar()
t_val.set('0')

output_dir = StringVar()
output_dir.set(r"D:\Manuscripts\2024_MotifQuest\MotifQuest_for_EG\output")

clustal_path = StringVar()
clustal_path.set(r"C:\Users\lawashburn\Downloads\clustal-omega-1.2.2\clustal-omega-1.2.2-win64\clustalo.exe")

def build_motif_db():

    from MotifQuest_code_v02 import start_building_motif_db
    start_building_a_motif_db()

def browse_files_fasta():
    filename = filedialog.askopenfilename(filetypes=[("FASTA Files",("*.fasta"))])
    fasta_path_input.set(filename)
    
def browse_files_clustalo():
    filename = filedialog.askopenfilename()
    clustal_path.set(filename)
    
def browse_files():
    filename = filedialog.askdirectory()
    output_dir.set(filename)

canvas = Canvas(
    window,
    bg = "#423C56",
    height = 750,
    width = 778,
    bd = 0,
    highlightthickness = 0,
    relief = "ridge"
)

canvas.place(x = 0, y = 0)
canvas.create_rectangle(
    11.0,
    101.0,
    748.0,
    517.0,
    fill="#D9D9D9",
    outline="")

canvas.create_text(
    17.0,
    123.0,
    anchor="nw",
    text="Import FASTA Database",
    fill="#000000",
    font=("Inter", 16 * -1)
)

canvas.create_text(
    5.0,
    0.0,
    anchor="nw",
    text="MotifQuest",
    fill="#FFFFFF",
    font=("Inter", 64 * -1)
)

entry_image_1 = PhotoImage(
    file=relative_to_assets("entry_1.png"))
entry_bg_1 = canvas.create_image(
    383.0,
    138.0,
    image=entry_image_1
)
entry_1 = Entry(
    bd=0,
    bg="#FFFFFF",
    highlightthickness=0,
    textvariable=fasta_path_input
)
entry_1.place(
    x=213.0,
    y=123.0,
    width=340.0,
    height=28.0
)

entry_image_2 = PhotoImage(
    file=relative_to_assets("entry_2.png"))
entry_bg_2 = canvas.create_image(
    383.0,
    224.0,
    image=entry_image_2
)
entry_2 = Entry(
    bd=0,
    bg="#FFFFFF",
    highlightthickness=0,
    textvariable=min_motif_len_in
)
entry_2.place(
    x=213.0,
    y=209.0,
    width=340.0,
    height=28.0
)

entry_image_3 = PhotoImage(
    file=relative_to_assets("entry_3.png"))
entry_bg_3 = canvas.create_image(
    383.0,
    274.0,
    image=entry_image_3
)
entry_3 = Entry(
    bd=0,
    bg="#FFFFFF",
    highlightthickness=0,
    textvariable=min_num_motif_inst
)
entry_3.place(
    x=213.0,
    y=259.0,
    width=200.0,
    height=28.0
)

entry_image_4 = PhotoImage(
    file=relative_to_assets("entry_4.png"))
entry_bg_4 = canvas.create_image(
    383.0,
    323.0,
    image=entry_image_4
)
entry_4 = Entry(
    bd=0,
    bg="#FFFFFF",
    highlightthickness=0,
    textvariable=part_motif_flank
)
entry_4.place(
    x=213.0,
    y=308.0,
    width=340.0,
    height=28.0
)

entry_image_5 = PhotoImage(
    file=relative_to_assets("entry_5.png"))
entry_bg_5 = canvas.create_image(
    383.0,
    370.0,
    image=entry_image_5
)
entry_5 = Entry(
    bd=0,
    bg="#FFFFFF",
    highlightthickness=0,
    textvariable = clustal_path
)
entry_5.place(
    x=213.0,
    y=400.0,
    width=340.0,
    height=28.0
)


button_image_1 = PhotoImage(
    file=relative_to_assets("motif_button_1.png"))
button_1 = Button(
    image=button_image_1,
    borderwidth=0,
    highlightthickness=0,
    command=lambda:browse_files_fasta(),
    relief="flat"
)
button_1.place(
    x=575.0,
    y=123.0,
    width=78.21710205078125,
    height=30.0
)

canvas.create_text(
    37.0,
    360.0,
    anchor="nw",
    text="Export Directory:",
    fill="#000000",
    font=("Inter", 16 * -1)
)

# entry_image_7 = PhotoImage(
#     file=relative_to_assets("entry_7.png"))
# entry_bg_7 = canvas.create_image(
#     383.0,
#     416.0,
#     image=entry_image_7
# )
entry_7 = Entry(
    bd=0,
    # bg="#FFFFFF",
    highlightthickness=0,
    textvariable=output_dir
)
entry_7.place(
    x=213.0,
    y=355.0,
    width=340.0,
    height=28.0
)

button_image_3 = PhotoImage(
    file=relative_to_assets("motif_button_3.png"))
button_3 = Button(
    image=button_image_3,
    borderwidth=0,
    highlightthickness=0,
    command=lambda: browse_files(),
    relief="flat"
)
button_3.place(
    x=575.0,
    y=355.0,
    width=78.21710205078125,
    height=30.0
)

canvas.create_text(
    33.0,
    173.0,
    anchor="nw",
    text="T-Value range",
    fill="#000000",
    font=("Inter", 16 * -1)
)

canvas.create_text(
    21.0,
    218.0,
    anchor="nw",
    text="Minimum Motif Length",
    fill="#000000",
    font=("Inter", 16 * -1)
)

canvas.create_text(
    15.0,
    260.0,
    anchor="nw",
    text="Minimum Number of\nMotif Instances",
    fill="#000000",
    font=("Inter", 16 * -1)
)

# entry_image_8 = PhotoImage(
#     file=relative_to_assets("entry_8.png"))
# entry_bg_8 = canvas.create_image(
#     249.0,
#     181.0,
#     image=entry_image_8
# )
entry_8 = Entry(
    bd=0,
    bg="#FFFFFF",
    highlightthickness=0,
    textvariable = t_val
)
entry_8.place(
    x=213.0,
    y=166.0,
    width=72.0,
    height=28.0
)

canvas.create_text(
    18.0,
    315.0,
    anchor="nw",
    text="Partial Motif Flank Size",
    fill="#000000",
    font=("Inter", 16 * -1)
)

canvas.create_text(
    18.0,
    411.0,
    anchor="nw",
    text="Aligned FASTA File Output",
    fill="#000000",
    font=("Inter", 16 * -1)
)

button_image_4 = PhotoImage(
    file=relative_to_assets("motif_button_4.png"))
button_4 = Button(
    image=button_image_4,
    borderwidth=0,
    highlightthickness=0,
    command= build_motif_db,
    relief="flat"
)
button_4.place(
    x=250.0,
    y=450.0,
    width=141.0,
    height=30.0
)

button_image_5 = PhotoImage(
    file=relative_to_assets("motif_button_5.png"))
button_5 = Button(
    image=button_image_5,
    borderwidth=0,
    highlightthickness=0,
    command = lambda: browse_files_clustalo(),
    relief="flat"
)
button_5.place(
    x=575.0,
    y=400.0,
    width=78.21710205078125,
    height=30.0
)
window.resizable(False, False)
window.mainloop()

