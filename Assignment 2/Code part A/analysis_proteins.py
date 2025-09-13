# -*- coding: utf-8 -*-
"""
Created on Sat Nov 12 16:53:31 2022

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

print("length of the reference sequence", len(ref_seq))
print("length of the given sequence", len(true_seq))

match = 1
mismatch = -2.5
opengap = -5
extendgap = -2
ps="match={0}, mismatch={1}, opengap={2}, extendgap={3}"
fps=ps.format(match, mismatch, opengap, extendgap)
scores=(match, mismatch, opengap, extendgap) 

cds = {} #allocate dictionary
#add keys to dictionay with as values lists of the two sliced sequences.
cds['spike protein'] = [ref_seq.seq[21563-1:25384], true_seq.seq[21000:26000]]
cds['orf3a'] = [ref_seq.seq[25393-1:26220], true_seq.seq[2500:27000]]
cds['envelope protein'] = [ref_seq.seq[26245-1:26472], true_seq.seq[25000:27000]]
cds['membrane glycoprotein'] = [ref_seq.seq[26523-1:27191], true_seq.seq[26000:28000]]
cds['orf6'] = [ref_seq.seq[27202-1:27387], true_seq.seq[2700:28000]]
cds['orf7a'] = [ref_seq.seq[27394-1:27759], true_seq.seq[27000:28000]]
cds['orf7b'] = [ref_seq.seq[27756-1:27887], true_seq.seq[27000:28000]]
cds['orf8'] = [ref_seq.seq[27894-1:28259], true_seq.seq[27000:29000]]
cds['nucleocapsid phosphoprotein'] = [ref_seq.seq[28274-1:29533], true_seq.seq[2900:30000]]
cds['orf10'] = [ref_seq.seq[29558-1:29674], true_seq.seq[29000:30000]]

protein_names = ['spike protein', 'orf3a', 'envelope protein', 'membrane glycoprotein', 'orf6', 'orf7a',
                  'orf7b', 'orf8', 'nucleocapsid phosphoprotein', 'orf10']

 
result_proteins = {}
for i in range(len(protein_names)): #loop over all the proteins
    protein_name = protein_names[i]
    #do semi-global alignment
    alignments = pairwise2.align.globalms(cds[protein_name][0], cds[protein_name][1], *scores,
                                          score_only=False,
                                          penalize_end_gaps=(False, False))
    nr_matches = 0
    for j in range(len(alignments[0][0])):
        if alignments[0][0][j] == alignments[0][1][j]: #check if character is the same
            nr_matches += 1 #if characters ar the same increase number of matches by 1
            
    #calculate the percentage of nucleotide with a mach
    perc_correct = (nr_matches/len(cds[protein_name][0]))*100 
    #add results to the dictionary
    result_proteins[protein_name] = [nr_matches, perc_correct]



