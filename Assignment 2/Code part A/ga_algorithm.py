# -*- coding: utf-8 -*-
"""
Created on Fri Nov 11 08:28:56 2022

@author: 20182460
"""

import numpy as np 

def globalAlignment(x, y, match, mismatch, opengap, extendgap):
    """
    gives a possible solution of global alignment where the allignment score is maximal. 

    Parameters
    ----------
    x: string with the first sequence
    y: string with the second sequence 
    match: score when there is a match 
    mismatch: score when there is a mismatch
    opengap: score when you open a gap (first gap after match or mismatch)
    extendgap: score when you extend a gap

    Returns
    -------
    final_score: the optimum global alignment score for the two sequences
    align_x: a possible optimal alignment for sequence x. 
    align_y: a possible optimal alignment for sequence y.

    """
    
    m=len(x)+1; 
    n=len(y)+1
    A= np.zeros((m, n))
    
    color_dict = {} #allocate dictionary where for each position in the score 
    #matrix a color and score is added
    for i in range(m): #loop over the number of number of characters is sequence x 
        #start with an open gap
        if i == 1:
            A[i,0] = opengap
            color_dict[(i,0)] = ('white', extendgap)
        #then extend gaps horizontally on first row
        else: 
            A[i,0] = A[i-1,0] + extendgap
            color_dict[(i,0)] =  ('white', extendgap)
    for j in range(n): #loop over the number of number of characters is sequence y
        #start with an open gap
        if j == 1:
            A[0,j] = opengap
            color_dict[(0,j)] = ('white', extendgap)
        #then extend gaps vertically on first colum
        else:
            A[0,j] = A[0,j-1] + extendgap
            color_dict[(0,j)] = ('white', extendgap)
    #make the very first entry zero
    A[0,0] = 0
    color_dict[(0,0)] = ('white', extendgap)
    
    for i in range(m-1):
        for j in range(n-1):
            if x[i] == y[j]:
                
                #find score coming from the tree possible directions
                left = A[i, j+1] + color_dict[(i,j+1)][1]
                diagonal = A[i, j] + match
                down =  A[i+1, j] + color_dict[(i+1,j)][1]
                #find maximum score and put this is new position
                max_score_match = max(left, diagonal, down)
                A[i+1, j+1] = max_score_match
                
                #add new position with color and score to the dictionary
                #when there was a match, the color becomes green and the score 
                #an opengap
                if max_score_match == diagonal and max_score_match != left and max_score_match != down:
                    color_dict[(i+1, j+1)] = ('green', opengap)
                else:
                    color_dict[(i+1, j+1)] = ('white', extendgap)
                   
            else:
                #find score coming from the tree possible directions
                left = A[i, j+1] + color_dict[(i,j+1)][1]
                diagonal = A[i, j] + mismatch
                down =  A[i+1, j] + color_dict[(i+1,j)][1]
                #find maximum score and put this is new position
                max_score_mis = max(left, diagonal, down)
                A[i+1, j+1] = max_score_mis
                #add new position with color and score to the dictionary
                #when there was a mismatch, the color becomes red and the score 
                #an opengap
                if max_score_mis == diagonal and max_score_mis != left and max_score_mis != down:
                    color_dict[(i+1, j+1)] = ('red', opengap)
                else:
                    color_dict[(i+1, j+1)] = ('white', extendgap)
                
                
    #tracing back scores by Needle-Wunsch principle
    align_x = ""
    align_y = ""
    
    i = len(x)
    j = len(y)
    
    while i > 0 or j > 0:
        print(i, j)
                 
        current_score = A[i][j]
            
        if x[i-1] == y[j-1]:
             similarity =  match
        else: 
             similarity =  mismatch
            
        # go to the left upper diagonal when there is a match or mismatch
        #and when the current position is green or red
        if i > 0 and j > 0 and current_score == A[i-1, j-1] + similarity \
        and (color_dict[(i,j)][0] == 'red' or color_dict[(i, j)][0] == 'green'):
            align_x = x[i - 1] + align_x
            align_y = y[j - 1] + align_y 
            i = i - 1
            j = j - 1
        
        # go up 
        elif i > 0 and (current_score == A[i-1, j] + extendgap or current_score == A[i-1, j] + opengap):
            align_x = x[i - 1] + align_x
            align_y = "-" + align_y
            i = i - 1
            
        # go to the left
        else:
            align_x = "-" + align_x
            align_y = y[j - 1] + align_y
            j = j - 1
        
    final_score = A[m-1, n-1]
    print('final score =', final_score)
    print(align_x)
    print(align_y)
        
    return final_score, align_x, align_y

x = 'ACCAATTACCAATTAAG'
y = 'AATGA'
globalAlignment(x, y, 1, -2.5, -5, -2)

print('ACCAATTACCAATTAAG')
print('----------AATGA--')

print('ACCAATTACCAATTAAG')
print('---AATGA---------')

print('ACCAATTACCAATTAAG')
print('A---ATGA---------')
