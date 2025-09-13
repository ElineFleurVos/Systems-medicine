# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 09:54:31 2022

@author: 20182460
"""


from Bio import Entrez
from Bio import SeqIO

def importSequence(filename):
    
    """
    imports sequence both in 'genbank' and 'fasta' format 

    Parameters
    ----------
    filename: string with the name of the sequence

    Returns
    -------
    seq_fasta: sequence in fasta format from the biopython module
    seq_gbk: sequence in genbank format from the biopython module 

    """

    Entrez.email = "e.f.vos@student.tue.nl"
    
    filename_gbk = filename + '.gbk'
    filename_fasta = filename + '.fa'
    
    #obtain sequence in genbank format and save it in the current folder
    net_handle_g = Entrez.efetch(db="nucleotide", id=filename, rettype="gb", retmode="text")
    out_handle_g = open(filename_gbk, 'w')
    out_handle_g.write(net_handle_g.read())
    out_handle_g.close()
    net_handle_g.close()
    #load sequence 
    seq_gbk = SeqIO.read(filename_gbk, 'genbank')
        
    #obtain sequence in fasta format and save it in the current folder
    net_handle_f = Entrez.efetch(db="nucleotide", id=filename, rettype="fasta", retmode="text")
    out_handle_f = open(filename_fasta, 'w')
    out_handle_f.write(net_handle_f.read())
    out_handle_f.close()
    net_handle_f.close()
    #load sequence 
    seq_fasta = SeqIO.read(filename_fasta, 'fasta')
    
    return seq_fasta, seq_gbk 
