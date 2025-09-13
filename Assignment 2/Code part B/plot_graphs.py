# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 09:12:57 2022

@author: 20182460
"""

from create_graphs import *
import cobra
import networkx as nx
import matplotlib.pyplot as plt

#%%
#load in the data 
#SET CORRECT PATH NAME
path_name = 'C:/Users/20182460/Desktop/Quartile 1/Systems medicine/Assignment 2/Assignment 2 8CM00 Eline Vos/'
m_aero1 = cobra.io.read_sbml_model(path_name + '/Bacteria/Aeromonas_jandaei_Riv2.xml')
m_aero2 = cobra.io.read_sbml_model(path_name + '/Bacteria/Aeromonas_hydrophila_subsp_hydrophila_ATCC_7966.xml')
m_clos = cobra.io.read_sbml_model(path_name + '/Bacteria/Clostridium_perfringens_ATCC_13124.xml')

#%%
def drawGraph(G, m):
    """
    Draws the graph G. A cobrypy metabolic model is also needed as input to assign the correct color per 
    node 

    Parameters
    ----------
    G : networkx graph
    m : Cobrapy metabolic model

    Returns
    -------
    None.

    """
    color_map = []
    for node in G:
        #make metabolites blue
        if node in m.metabolites:
            color_map.append('royalblue')
        #make reactions red
        else: 
            color_map.append('red')   
            
    plt.figure()
    options = {
        'node_size': 50,
        'width': 3 }
    nx.draw(G, **options, node_color = color_map)
    plt.show()
    
def plotGraphs(m):
    """
    plots the graphs and largest subgraph for graph made of m using drawGraph.

    Parameters
    ----------
    m : Cobrapy metabolic model

    Returns
    -------
    None.

    """
    Gm, size_Gm = createMetaboliteGraph(m)
    Gmr, size_Gmr = createMetaReacGraph(m)
    
    drawGraph(Gm, m)
    drawGraph(Gmr, m)
    
    Sm = [Gm.subgraph(c).copy() for c in nx.connected_components(Gm)]
    Smr = [Gmr.subgraph(c).copy() for c in nx.connected_components(Gmr)]
    
    drawGraph(Sm[0], m)
    drawGraph(Smr[0], m)
    
    return None

#%%

plotGraphs(m_aero1)
plotGraphs(m_aero2)
plotGraphs(m_clos)
    