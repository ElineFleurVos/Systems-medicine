# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 00:00:10 2022

@author: 20182460
"""

import cobra
import pandas as pd
import matplotlib.pyplot as plt
from cobra.util.solver import linear_reaction_coefficients
from shortest_path import *
from create_graphs import *

#%%

def makeEmptyMedium(m):
    """
    changes the medium of the cobrapymodel to 'empty' by setting all flux bounds to zero.

    Parameters
    ----------
    m : Cobrypy metabolic model

    Returns
    -------
    m : Cobrapy metabolic model where the medium is 'empty'
    
    """
    medium = m.medium
    medium_names = list(medium.keys())
    for i in range(len(medium_names)): #loop over all medium metabolites
        medium_name = medium_names[i]
        medium[medium_name] = 0.0 #set flux to zero
    
    #set new medium
    m.medium = medium
    
    return m

def makeDMEM(m):
    """
    changes the medium of the cobrapymodel to a DMEM medium by setting most of the flux bound to zero except 
    for EX_gly_L(e), EX_glc_D(e), ES_na1(e)
    
    Parameters
    ----------
    m : Cobrapy metabolic model

    Returns
    -------
    m : Cobrapy metabolic model with a DMEM medium

    """
    medium_DMEM = ["EX_glu_L(e)", "EX_glc_D(e)", "EX_na1(e)"]
    # L_glutamine = m.metabolites.get_by_id("glu_L_e")
    # glucose = m.metabolites.get_by_id("glc_D_e")
    # sodium = m.metabolites.get_by_id("na1_e")
            
    medium = m.medium
    medium_names = list(medium.keys()) #get names of metabolites in medium
    for i in range(len(medium_names)): #loop over all metabolites in medium
        medium_name = medium_names[i]
        if medium_name not in medium_DMEM: #check if one of the metabolites in DMEM medium
            medium[medium_name] = 0.0
    
    #set new medium
    m.medium = medium
    
    return m
     
def degreeMediumMetabolites(G, m):
    """
    makes a dictionary with names of the metabolites in the medium as key and as value the degree they 
    have in the metabolite graph (degree of the variant in the compartment and not the extracellular space)

    Parameters
    ----------
    G : NetworkX metabolite graph
    m : Cobrapy metabolic model from which the graph is made

    Returns
    -------
    medium_metabolite_degree : dictionary with as keys the names of the metabolites in the medium of m and 
                               as values the degree in the metabolite graph.
    """
    
    mediumG = m.medium #get the exchange reactions with their fluxes.
    medium_metabolites = []
    #get all metabolites in the medium 
    for j in range(len(m.metabolites)):
        metabolite = m.metabolites[j]
        if metabolite.compartment == 'e':
            medium_metabolites.append(metabolite.id) 
    
    #this part is added because for Aero2, the medium_metabolites list above has more metabolites 
    #then wen you use m.medium directly from the library. The code did not work in that case.
    if len(m.metabolites) == 1280:
        medium_metabolites.remove('tmao_e')
        medium_metabolites.remove('tet_e')
        medium_metabolites.remove('no3_e')
        medium_metabolites.remove('no2_e')
    
    #change the names to get the metabolites in the bacterium itself, because those degrees are more useful
    compartment_meta = []
    for i in medium_metabolites:
        mm = i[:len(i)-1]
        mc = mm + 'c'
        compartment_meta.append(mc)
              
    medium_metabolite_degree = {}
    list_medium_keys = list(mediumG.keys())
    for i in range(len(compartment_meta)):
        meta = compartment_meta[i]
        degree = G.degree(meta) #get the degree for the metabolite 
        if isinstance(degree, int):
            medium_metabolite_degree[list_medium_keys[i]] = degree #add name and degree to dictionary
          
    return medium_metabolite_degree

def addMetabolites(m, medium_degrees):
    """
    This function systematically adds metabolites to the medium by adding the metabolites with the
    highest degree one by one until the bacterium is able to grow in the medium

    Parameters
    ----------
    m : Cobrapy metabolic model
    medium_degrees : dictionary with as keys the names of the metabolites in the medium of m and 
                     as values the degree in the metabolite graph.

    Returns
    -------
    added : list with the metabolites that were added to the medium.
    degree_added : list with the degrees of the metabolites that were added to the medium.
    final_score : the first objective value score found that is higher dan 0.  

    """

    score = 0
    added = []
    degree_added = []
    medium = m.medium
    while score == 0: #go on until objective value is bigger than zero
        medium_names = list(medium_degrees.keys())
        max_key = max(medium_degrees, key=medium_degrees.get) #get the name of the metabolite with the highest degree
        max_value = medium_degrees[max_key] #get highest degree
        for i in medium_names:
            if i == max_key:
                medium[i] = 1000 #set upper bound back to 0 for metabolite with highest degree
                added.append(i)
                degree_added.append(max_value)
                del medium_degrees[i] #delete the metabolite, so that the next max degree can be found
        
        #set new medium and calculate the new objective value
        m.medium = medium
        s = m.optimize()
        score = s.objective_value
    
    #get final objective value score
    s_final = m.optimize()
    final_score = s_final.objective_value
    
    return added, degree_added, final_score