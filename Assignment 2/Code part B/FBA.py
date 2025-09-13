# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 09:30:37 2022

@author: 20182460
"""

import cobra
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from cobra.util.solver import linear_reaction_coefficients
from create_graphs import createWeightedGraph
import matplotlib.pyplot as plt
import math
import operator
from matplotlib.ticker import ScalarFormatter
from shortest_path import *

#%%
#load in the data 
#SET CORRECT PATH NAME
path_name = 'C:/Users/20182460/Desktop/Quartile 1/Systems medicine/Assignment 2/Assignment 2 8CM00 Eline Vos/'
m_aero1 = cobra.io.read_sbml_model(path_name + '/Bacteria/Aeromonas_jandaei_Riv2.xml')
m_aero2 = cobra.io.read_sbml_model(path_name + '/Bacteria/Aeromonas_hydrophila_subsp_hydrophila_ATCC_7966.xml')
m_clos = cobra.io.read_sbml_model(path_name + '/Bacteria/Clostridium_perfringens_ATCC_13124.xml')
    
#%%

def get_biomass_reaction(m):
    """
    Return the biomass reaction of a model

    Parameters
    ----------
    m : Cobrapy metabolic model

    Returns
    -------
    rec : biomass reaction of the model

    """
    objective_str = str(list(m.objective.variables)[0])
    for rec in m.reactions:
        if rec.id in objective_str:
            return rec

#get the biomass reaction of the three bacteria
biomass_aero1 = get_biomass_reaction(m_aero1)
biomass_aero2 = get_biomass_reaction(m_aero2)
biomass_closs = get_biomass_reaction(m_clos)

#get the objective value for each of the bacteria
biomass_aero1 = m_aero1.reactions.get_by_id("biomass423")
m_aero1.objective = biomass_aero1
s1 = m_aero1.optimize()
score1 = s1.objective_value

biomass_aero1 = m_aero2.reactions.get_by_id("biomass444")
m_aero2.objective = biomass_aero2
s2 = m_aero2.optimize()
score2 = s2.objective_value

biomass_clos = m_clos.reactions.get_by_id("biomass012")
m_clos.objective = biomass_clos
s_clos = m_clos.optimize()
score_clos = s_clos.objective_value

#get all shortest paths of the metabolites in the medium to the biomass
min_value1, max_value1, list_lengths1, list_freq1 = shortestPath(m_aero1)
min_value2, max_value2, list_lengths2, list_freq2 = shortestPath(m_aero2)
min_value3, max_value3, list_lengths3, list_freq3 = shortestPath(m_clos)

#plot the shortest path lengths with the highest frequency
ax = plotMostFrequent(m_aero1, m_aero2, m_clos)
add_value_labels(ax)   
    
