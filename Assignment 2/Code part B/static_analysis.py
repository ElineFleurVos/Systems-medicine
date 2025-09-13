# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 14:45:32 2022

@author: 20182460
"""

import cobra
import networkx as nx
from create_graphs import *
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
import pandas as pd
import operator

#%%
#load in the data 
#SET CORRECT PATH NAME
path_name = 'C:/Users/20182460/Desktop/Quartile 1/Systems medicine/Assignment 2/Assignment 2 8CM00 Eline Vos/'
m_aero1 = cobra.io.read_sbml_model(path_name + '/Bacteria/Aeromonas_jandaei_Riv2.xml')
m_aero2 = cobra.io.read_sbml_model(path_name + '/Bacteria/Aeromonas_hydrophila_subsp_hydrophila_ATCC_7966.xml')
m_clos = cobra.io.read_sbml_model(path_name + '/Bacteria/Clostridium_perfringens_ATCC_13124.xml')

#%%

#create all the graphs 
Gm_aero1, size_Gm_aero1 = createMetaboliteGraph(m_aero1)
Gmr_aero1, size_Gmr_aero1 = createMetaReacGraph((m_aero1))
Gm_aero2, size_Gm_aero2 = createMetaboliteGraph(m_aero2)
Gmr_aero2, size_Gmr_aero2 = createMetaReacGraph((m_aero2))
Gm_clos, size_Gm_clos = createMetaboliteGraph(m_clos)
Gmr_clos, size_Gmr_clos = createMetaReacGraph((m_clos))

#%% make degree distribution plots 

#sort degrees of nodes in metabolite graphs
degree_sequence1 = sorted((d for n, d in Gm_aero1.degree()), reverse=True)
degree_sequence2 = sorted((d for n, d in Gm_aero2.degree()), reverse=True)
degree_sequence_clos = sorted((d for n, d in Gm_clos.degree()), reverse=True)

fig, ax = plt.subplots(1,3)
fig.set_size_inches(20,4)
ax[0].axis([1, 2000, 1, 2000])
#plot loglog vs rank of degrees
ax[0].loglog(degree_sequence1, marker="o", base = 2)
ax[1].loglog(degree_sequence2, marker="o", base = 2)
ax[2].loglog(degree_sequence_clos, marker="o", base = 2)
for i in range(3):
    ax[i].set_xlabel('rank', fontsize =  16)
    ax[i].set_ylabel('metabolite degree',fontsize =  16)
    ax[i].set_xticks((1, 5, 10, 20, 50, 200, 400, 1000))
    ax[i].set_yticks((1, 5, 10, 20, 50, 200, 400, 1200))
    for axis in [ax[i].xaxis, ax[i].yaxis]:
        axis.set_major_formatter(ScalarFormatter())
    for tick in ax[i].xaxis.get_major_ticks():
                tick.label.set_fontsize(14) 
    for tick in ax[i].yaxis.get_major_ticks():
                tick.label.set_fontsize(14) 

#%% find hubs 

def findHubs(G):
    """
    Find the number of hubs in metabolite graph G. 
    The definition of the hubs can be found in the report (Equation 1)

    Parameters
    ----------
    G : NetworkX metabolite graph 

    Returns
    -------
    degree_results : dictionary with as keys the names of the metabolites and as value the degree of the 
                     metabolite
    nr_hubs : integer of the amount of hubs

    """
    all_degrees = []
    degree_results = {}
    #get the degree for all nodes and put them in a list and dictionary
    for i in G.nodes:
        degree_results[i] = G.degree(i)
        all_degrees.append(G.degree(i))
        
    #calculate sigma
    avg_degree = np.mean(all_degrees)
    N = len(Gm_aero1)
    p = avg_degree/(N-1)
    sigma2 = (N-1)*p*(1-p)
    sigma = np.sqrt(sigma2)

    #assign metabolite to hubs list if degree larger than mean degree plus three times sigma
    hubs = []
    for node in degree_results:
        degree = degree_results[node]
        if degree >= avg_degree + 3*sigma:
            hubs.append(node)
    nr_hubs = len(hubs)
          
    #plot 10 metabolites with highest degree
    high_keys = []
    high_degrees = []
                
    #get 10 metabolites with highest degree in put names in high_keys and degrees in high_degrees
    highest_degrees = dict(sorted(degree_results.items(), key=operator.itemgetter(1), reverse=True)[:10])
    for i in highest_degrees:
        high_keys.append(i)
        high_degrees.append(highest_degrees[i])
        
    #make bar plot
    plt.figure()
    plt.bar(high_keys, high_degrees)
    plt.xlabel('metabolite')
    plt.ylabel('degree')
    
    return degree_results, nr_hubs

hubs_aero1, nr_hubs_aero1 = findHubs(Gm_aero1)
hubs_aero2, nr_hubs_aero2 = findHubs(Gm_aero2)
hubs_clos, nr_hubs_clos = findHubs(Gm_clos)


