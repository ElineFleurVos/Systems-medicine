# -*- coding: utf-8 -*-
"""
Created on Sun Nov 13 10:44:53 2022

@author: 20182460
"""


from load_sequence import importSequence
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from translation_function import translate

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

spike_protein = [ref_seq.seq[21563-1:25384], true_seq.seq[21000:26000]]

#do global alignment for spike protein 
alignments = pairwise2.align.globalms(spike_protein[0], spike_protein[1], *scores,
                                      score_only=False,
                                      penalize_end_gaps=(False, False))

#find the start positon of the coding sequence on the provided sequence
for i in range(len(alignments[0][0])):
    if alignments[0][0][i] != '-':
        start = i
        break
#find the end position of the coding sequence on the provided sequence
for i in range(len(alignments[0][0])):
    irev = len(alignments[0][0])-i-1
    if alignments[0][0][irev] != '-':
       end = irev
       break

#make a string with the reference RNA of the spike protein 
spike_ref_RNA = ''
for i in range(len(alignments[0][0])):
    if alignments[0][0][i] != '-':  #remove gaps
        spike_ref_RNA += alignments[0][0][i]

#make a string with the RNA of the spike protein on the provided sequence.
spike_RNA = alignments[0][1][start:end+1] #slice the correct part of the provided sequence
spike_true_RNA = ""
for i in range(len(spike_RNA)):
    if spike_RNA[i] != '-':   #remove gaps
        spike_true_RNA += spike_RNA[i]

#translate the RNA sequences to amino acid sequences
spike_ref_amino = translate(spike_ref_RNA)
spike_true_amino = translate(spike_true_RNA)

#do global alignment to align perfectly. Than you can see where the deletions are on the 
#provided amino acid sequence
match_s = 1
mismatch_s = -1
opengap_s = -2
extendgap_s = -2
ps="match={0}, mismatch={1}, opengap={2}, extendgap={3}"
fps_spike=ps.format(match, mismatch, opengap, extendgap)
scores_spike =(match_s, mismatch_s, opengap_s, extendgap_s)

alignments_spike = pairwise2.align.localms(spike_ref_amino, spike_true_amino, *scores_spike,
                                           score_only=False)

#%% look at specific positions based on covariants.org and put results in dictionary

positions_deduced = [19,24,25,26,69, 144, 157, 452]
positions = [19,24,25,26,27,142,213,339,371,373,375,376,405,408,417,440,452,477,478,484,493,498,501,505,614,655,679,681,704,764,796,954,969]
spike_result = {}
for i in range(len(positions)):
    position = positions[i]
    ref_amino = alignments_spike[0][0][position-1]
    true_amino = alignments_spike[0][1][position-1]
    spike_result[position] = (ref_amino, true_amino)
        



