# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 15:30:07 2022

@author: 20182460
"""

from load_sequence import importSequence
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

#%% load in all data 
ref_seq = 'NC_045512'
ref_fasta, ref_genbank = importSequence(ref_seq)
ref_seq = ref_fasta
true_seq = SeqIO.read('vos.fa', 'fasta')

#%%

match = 1
mismatch = -2.5
opengap = -5
extendgap = -2
ps="match={0}, mismatch={1}, opengap={2}, extendgap={3}"
fps=ps.format(match, mismatch, opengap, extendgap)
scores=(match, mismatch, opengap, extendgap) 

#do global alignment for membrane glycoprotein and look manually what is wrong
membrane_glycoprotein = [ref_seq.seq[26523-1:27191], true_seq.seq[26000:28000]]

alignments = pairwise2.align.globalms(membrane_glycoprotein[0], membrane_glycoprotein[1], *scores,
                                      score_only=False,
                                      penalize_end_gaps=(False, False))

#count the number of N's in the provided sequence
amount_N = 0
for i in range(len(alignments[0][1])):
    if alignments[0][1][i] == 'N':
        amount_N += 1
    

