# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 10:37:07 2022

@author: 20182460
"""

#FBA analyis starting with DMEM as medium 

import cobra
import pandas as pd
import matplotlib.pyplot as plt
from cobra.util.solver import linear_reaction_coefficients
from shortest_path import *
from create_graphs import *
from change_media import *

#%% START WITH DMEM MEDIUM
#SET CORRECT PATH NAME
path_name = 'C:/Users/20182460/Desktop/Quartile 1/Systems medicine/Assignment 2/Assignment 2 8CM00 Eline Vos/'
m_aero1 = cobra.io.read_sbml_model(path_name + '/Bacteria/Aeromonas_jandaei_Riv2.xml')
m_aero2 = cobra.io.read_sbml_model(path_name + '/Bacteria/Aeromonas_hydrophila_subsp_hydrophila_ATCC_7966.xml')
m_clos = cobra.io.read_sbml_model(path_name + '/Bacteria/Clostridium_perfringens_ATCC_13124.xml')

#%% 
#THIS ONLY WORK WHEN YOU RUN IN THIS ORDER. Function degree_medium_metabolite only works when function
#makeDMEM is not run yet. You can run it just once or you have to reset the variables.

#Run all the results when adding metabolites to the DMEM medium

Gm1, size = createMetaboliteGraph(m_aero1)
medium_degree1d = degreeMediumMetabolites(Gm1, m_aero1)
m_aero1_DMEM = makeDMEM(m_aero1)
added_dmem1, degree_added_dmem1, score_dmem1 = addMetabolites(m_aero1_DMEM, medium_degree1d)
min_value1d, max_value1d, list_lengths1d, list_freq1 = shortestPath(m_aero1_DMEM)

Gm2, size = createMetaboliteGraph(m_aero2)
medium_degree2d = degreeMediumMetabolites(Gm2, m_aero2)
m_aero2_DMEM = makeDMEM(m_aero2)
added_dmem2, degree_added_dmem2, score_dmem2 = addMetabolites(m_aero2_DMEM, medium_degree2d)
min_value2d, max_value2d, list_lengths2d, list_freq2d = shortestPath(m_aero2_DMEM)

Gm3, size = createMetaboliteGraph(m_clos)
medium_degree3d = degreeMediumMetabolites(Gm3, m_clos)
m_clos_DMEM = makeDMEM(m_clos)
added_dmem3, degree_added_dmem3, score_dmem3 = addMetabolites(m_clos_DMEM, medium_degree3d)
min_value3d, max_value3d, list_lengths3d, list_freq3d = shortestPath(m_clos_DMEM)

#plot most frequent path lengths
ax = plotMostFrequent(m_aero1, m_aero2, m_clos)
add_value_labels(ax)  

#%% START WITH EMPTY MEDIUM
#load in the data again because the code before chagned the bacteria
m_aero1 = cobra.io.read_sbml_model(path_name + '/Bacteria/Aeromonas_jandaei_Riv2.xml')
m_aero2 = cobra.io.read_sbml_model(path_name + '/Bacteria/Aeromonas_hydrophila_subsp_hydrophila_ATCC_7966.xml')
m_clos = cobra.io.read_sbml_model(path_name + '/Bacteria/Clostridium_perfringens_ATCC_13124.xml')

#%%
#THIS ONLY WORK WHEN YOU RUN IN THIS ORDER. Function degree_medium_metabolite only works when function
#makeEmptuMedium is not run yet. You can run it just once or you have to reset the variables.

#Run all the results when adding metabolites to the empty medium
Gm1, size = createMetaboliteGraph(m_aero1)
medium_degree1e = degreeMediumMetabolites(Gm1, m_aero1)
m_aero1_empty = makeEmptyMedium(m_aero1)
added_empty1, degree_added_empty1, score_empty1 = addMetabolites(m_aero1_empty, medium_degree1e)
min_value1e, max_value1e, list_lengths1e, list_freq1e = shortestPath(m_aero1_empty)

Gm2, size = createMetaboliteGraph(m_aero2)
medium_degree2e = degreeMediumMetabolites(Gm2, m_aero2)
m_aero2_empty = makeEmptyMedium(m_aero2)
added_empty2, degree_added_empty2, score_empty2 = addMetabolites(m_aero2_empty, medium_degree2e)
min_value2e, max_value2e, list_lengths2e, list_freq2e = shortestPath(m_aero2_empty)

Gm3, size = createMetaboliteGraph(m_clos)
medium_degree3e = degreeMediumMetabolites(Gm3, m_clos)
m_clos_empty = makeEmptyMedium(m_clos)
added_empty3e, degree_added_empty3, score_empty3 = addMetabolites(m_clos_empty, medium_degree3e)
min_value3e, max_value3e, list_lengths3e, list_freq3e = shortestPath(m_clos_empty)

#plot most frequent path lengths
ax = plotMostFrequent(m_aero1_empty, m_aero2_empty, m_clos_empty)
add_value_labels(ax)   
    

        
        
        
        
        