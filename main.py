#GUI
import tkinter as tk
#BIOINFORMATICS
from Bio import Entrez, SeqIO
#OTHERS
from collections import Counter
import ctypes
import os
import sys

base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))

dll_path = os.path.join(base_path, "compute.dll")

#variables
protein_aacids = ""
codon_usage_glb = ""

def calculate_codon_usage(sequences):
    """Calculeaza valorile pentru optimizarea codonilor."""
    if sequences is None: return None
    codon_counts = Counter()
    
    for seq in sequences:
        for i in range(0, len(seq) - 2, 3):  # Read codons (triplets)
            codon = seq[i:i+3].upper()  # Ensure uppercase codons
            codon_counts[codon] += 1
    
    total_codons = sum(codon_counts.values()) if codon_counts else 1
    codon_usage = {str(codon): round(count / total_codons * 100, 2) for codon, count in codon_counts.items()}
    
    return codon_usage

def fetch_coding_sequences(organism_name):
    if not organism_name: return None
    """Coding seq"""
    search_query = f"{organism_name}[Organism] AND CDS[Feature]"

    Entrez.email = "raresc589@gmail.com"
    
    search_handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=1)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    
    if not search_results["IdList"]:
        print(f"No coding sequences found for {organism_name}.")
        return None
    
    seq_id = search_results["IdList"][0]
    
    fetch_handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
    seq_record = SeqIO.read(fetch_handle, "genbank")
    fetch_handle.close()
    
    coding_sequences = []
    
    for feature in seq_record.features:
        if feature.type == "CDS":
            if "translation" in feature.qualifiers:
                coding_sequences.append(feature.extract(seq_record.seq))
    
    return coding_sequences

def get_aminoacid_sequence(protein_name, codon_usage, nb):
    global protein_aacids
    number=0

    try:
        number=int(nb)
    except:
        number=15

    Entrez.email = "raresc589@gmail.com"
    
    handle = Entrez.esearch(db="protein", term=protein_name, retmax=1)
    record = Entrez.read(handle)
    handle.close()

    
    if not record["IdList"] or protein_name.isupper():
        if (not protein_name.isupper()):
            set_text(f"Protein {protein_name} not found. Please enter protein sequence manually.\n")
        protein_aacids=protein_name

    else:

        protein_id = record["IdList"][0]

        handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
        seq_record = SeqIO.read(handle, "fasta")
        handle.close()

        set_text("Fetched protein sequence.\n")
        protein_aacids=str(seq_record.seq)

    l6.configure(text="Total number of nucleotide bases: "+str(len(protein_aacids)))

    if codon_usage:
        codon_usage_glb = codon_usage 
        set_text(f"Protein sequence: {protein_aacids}\n\nCodon usage:\n{codon_usage}\n")

        #codon_usage - hashmap
        #protein_aacids - string

        lib = ctypes.CDLL(dll_path)

        lib.get_combinations.argtypes = (ctypes.c_int,ctypes.c_char_p)
        lib.get_combinations.restype = ctypes.c_double
        protein_bytes = ctypes.create_string_buffer(protein_aacids.encode('utf-8'))

        l3.config(text="Total combinations: 10^"+str(lib.get_combinations(0, protein_bytes)))
        l4.config(text="Total optimised combinations: 10^"+str(lib.get_combinations(1, protein_bytes)))

        keys = (ctypes.c_char_p * len(codon_usage.keys()))(*[key.encode('utf-8') for key in codon_usage.keys()])
        values = (ctypes.c_float * len(codon_usage.values()))(*codon_usage.values())
        lib.compute.argtypes = (ctypes.POINTER(ctypes.c_char_p),ctypes.POINTER(ctypes.c_float),ctypes.c_char_p,ctypes.c_int)
        lib.compute.restype = None
        lib.compute(keys, values, protein_bytes, number)

    else:
        set_text("No codon usage data available.\n")
        set_text(f"Protein sequence: {protein_aacids}\n")

        lib = ctypes.CDLL(dll_path)
        
        lib.get_combinations.argtypes = (ctypes.c_int,ctypes.c_char_p)
        lib.get_combinations.restype = ctypes.c_double

        protein_bytes = ctypes.create_string_buffer(protein_aacids.encode('utf-8'))

        l3.config(text="Total combinations: 10^"+str(lib.get_combinations(0,protein_bytes)))
        l4.config(text="Total optimised combinations: --")

        lib.compute.argtypes = (ctypes.POINTER(ctypes.c_char_p),ctypes.POINTER(ctypes.c_float),ctypes.c_char_p,ctypes.c_int)
        lib.compute.restype = ctypes.c_char_p
        lib.compute(ctypes.POINTER(ctypes.c_char_p)(), ctypes.POINTER(ctypes.c_float)(), protein_bytes, number)





#gui
window = tk.Tk()
window.title("ADN Sequence Finder")
window.resizable(0,0)
window.geometry("500x500")

def open_popup(type):
    popup = tk.Toplevel(window)
    if type == "help":
        popup.title("Help")
        popup.geometry("300x200")
        text=open("./data/help.txt","r").read()
        l = tk.Label(popup, text=text)
        l.pack()
    elif type == "about":
        popup.title("About")
        popup.geometry("300x200")
        text=open("./data/about.txt","r").read()
        l = tk.Label(popup, text=text)
        l.pack()
    else:
        popup.title("Contact")
        popup.geometry("300x200")
        text=open("./data/contact.txt","r").read()
        l = tk.Label(popup, text=text)
        l.pack()

menu=tk.Menu(window)
window.config(menu=menu)
helpmenu=tk.Menu(menu,tearoff=0)
menu.add_cascade(label="Help",menu=helpmenu)
helpmenu.add_command(label="Help",command=lambda: open_popup("help"))
helpmenu.add_command(label="About",command=lambda: open_popup("about"))
helpmenu.add_command(label="Contact",command=lambda: open_popup("contact"))


l=tk.Label(window,text="Enter the protein name (or amino acid sequence in uppercase):")
l.pack()
e=tk.Entry(window)
e.pack()
l2=tk.Label(window,text="Enter the organism name:")
l2.pack()
e2=tk.Entry(window)
e2.pack()
b=tk.Button(window,text="Find",command=lambda: get_aminoacid_sequence(e.get(),calculate_codon_usage(fetch_coding_sequences(e2.get())),e3.get()))
b.pack()

l3=tk.Label(window,text="Total combinations: --")
l3.pack()
l4=tk.Label(window,text="Total optimised combinations: --")
l4.pack()
l6=tk.Label(window,text="Number of nucleotide bases: --")
l6.pack()
l5=tk.Label(window,text="Enter number of combinations to output:")
l5.pack()
e3=tk.Entry(window)
e3.pack()


frame = tk.Frame(window)
frame.pack(padx=10, pady=10, fill="both", expand=False)

text_widget = tk.Text(frame, wrap="word", height=10, width=75, font=("Times new roman", 10), state="disabled")
text_widget.pack(side="left", fill="both", expand=False)

scrollbar = tk.Scrollbar(frame, command=text_widget.yview)
scrollbar.pack(side="left", fill="y")

text_widget.config(yscrollcommand=scrollbar.set)

def set_text(content,clear=False,position="1.0"):
    text_widget.config(state="normal")
    if clear: text_widget.delete(position, tk.END) 
    text_widget.insert("1.0", content)
    text_widget.config(state="disabled")

window.mainloop()