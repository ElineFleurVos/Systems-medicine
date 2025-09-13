# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 09:54:46 2022

@author: 20182460
"""

#This file checks the difference between calculating the similarity score as a sum of two chunked 
#sequences compared to doing the semi global alignment in one go. Also,the similarity score and total 
#amount of matches is calculated when using the maximal chunk size of 5000. 

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
print("length of the reference sequence", len(ref_seq))
print("length of the given sequence", len(true_seq))

#define score function
match = 1
mismatch = -2.5
opengap = -5
extendgap = -2
ps="match={0}, mismatch={1}, opengap={2}, extendgap={3}"
fps=ps.format(match, mismatch, opengap, extendgap)
scores=(match, mismatch, opengap, extendgap) 

#semi global alignment of first 2000 nucleotides
ref_seq_short = ref_seq.seq[:2000]
true_seq_short = true_seq.seq[:2000]

scores=(match, mismatch, opengap, extendgap) 
alignments = pairwise2.align.globalms(ref_seq_short, true_seq_short, *scores,
                                      score_only=False,
                                      penalize_end_gaps=(False, False))

#semi global alignment of nucleotides 2001 until 4000
ref_seq_short2 = ref_seq.seq[2001:4000]
true_seq_short2 = true_seq.seq[2001:4000]

scores=(match, mismatch, opengap, extendgap) 
alignments2 = pairwise2.align.globalms(ref_seq_short2, true_seq_short2, *scores,
                                      score_only=False,
                                      penalize_end_gaps=(False, False))

#semi global alignment of first 4000 nucleotides
ref_seq_combine = ref_seq.seq[:4000]
true_seq_combine  = true_seq.seq[:4000]

scores=(match, mismatch, opengap, extendgap) 
alignments_combine  = pairwise2.align.globalms(ref_seq_combine, true_seq_combine, *scores,
                                      score_only=False,
                                      penalize_end_gaps=(False, False))

#print the chunked results and the results in one go
#length alignments
print("in chunks", len(alignments[0][0]) + len(alignments2[0][0]))
print("in one go", len(alignments_combine[0][0]))
#score alignments
print("in chunks", alignments[0][2] + alignments2[0][2])
print("in one go", alignments_combine[0][2])

#%% calculate total amount of matches when using chunks

n_start = 0
n_end = 5000

total_score = 0
total_matches = 0
nr_alignments = []
for i in range(6):
    
    ref_seq_short = ref_seq.seq[n_start:n_end]
    true_seq_short = true_seq.seq[n_start:n_end]

    scores=(match, mismatch, opengap, extendgap) 
    alignments = pairwise2.align.globalms(ref_seq_short, true_seq_short, *scores,
                                          score_only=False,
                                          penalize_end_gaps=(False, False))
    nr_alignments.append(len(alignments))
    total_score += alignments[0][2]
    n_start += 5000
    n_end += 5000
    
    for j in range(len(alignments[0][0])):
        if alignments[0][0][j] == alignments[0][1][j]:
            total_matches += 1
            
print("similarity score =", total_score)
print("number of matches =", total_matches)
    
    
    
    
    
    








