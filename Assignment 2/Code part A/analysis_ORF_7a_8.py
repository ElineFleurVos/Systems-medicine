# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 16:42:25 2022

@author: 20182460
"""

from load_sequence import importSequence
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from translation_function import translate

#this file checks if the mutations on ORF7a and ORF8 of the Delta 21A variant are 
#really present on the reference sequence NC_045512

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

#%% do analysis for ORF7a protein 

orf7a = [ref_seq.seq[27394-1:27759], true_seq.seq[27000:28000]]

alignments_7a = pairwise2.align.globalms(orf7a[0], orf7a[1], *scores,
                                      score_only=False,
                                      penalize_end_gaps=(False, False))

for i in range(len(alignments_7a[0][0])):
    if alignments_7a[0][0][i] != '-':
        start = i
        break
for i in range(len(alignments_7a[0][0])):
    irev = len(alignments_7a[0][0])-i-1
    if alignments_7a[0][0][irev] != '-':
       end = irev
       break

ref_RNA_7a = ''
for i in range(len(alignments_7a[0][0])):
    if alignments_7a[0][0][i] != '-':
        ref_RNA_7a += alignments_7a[0][0][i]
                
RNA_7a = alignments_7a[0][1][start:end+1]
true_RNA_7a = ""
for i in range(len(RNA_7a)):
    if RNA_7a[i] != '-':
        true_RNA_7a += RNA_7a[i]
        
ref_amino_7a = translate(ref_RNA_7a)
true_amino_7a = translate(true_RNA_7a)
        
#check if amino sequences are the same.         
nr_matches_7b = 0
for i in range(len(ref_amino_7a)):
    if ref_amino_7a[i] == true_amino_7a[i]:
        nr_matches_7b += 1

print(ref_amino_7a[81])
print(ref_amino_7a[119])

#%% do analyis for ORF8 protein 

orf8 = [ref_seq.seq[27894-1:28259], true_seq.seq[27000:29000]]

alignments_8 = pairwise2.align.globalms(orf8[0], orf8[1], *scores,
                                      score_only=False,
                                      penalize_end_gaps=(False, False))


for i in range(len(alignments_8[0][0])):
    if alignments_8[0][0][i] != '-':
        start = i
        break
for i in range(len(alignments_8[0][0])):
    irev = len(alignments_8[0][0])-i-1
    if alignments_8[0][0][irev] != '-':
       end = irev
       break

ref_RNA_8 = ''
for i in range(len(alignments_8[0][0])):
    if alignments_8[0][0][i] != '-':
        ref_RNA_8 += alignments_8[0][0][i]
                
RNA_8 = alignments_8[0][1][start:end+1]
true_RNA_8 = ""
for i in range(len(RNA_8)):
    if RNA_8[i] != '-':
        true_RNA_8 += RNA_8[i]
        
ref_amino_8 = translate(ref_RNA_8)
true_amino_8 = translate(true_RNA_8)
        
#check if amino sequences are the same.         
nr_matches_8 = 0
for i in range(len(ref_amino_8)):
    if ref_amino_8[i] == true_amino_8[i]:
        nr_matches_8 += 1

print(ref_amino_8[118])
print(ref_amino_8[119])
