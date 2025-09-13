# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 08:26:07 2022

@author: 20182460
"""

import networkx as nx

#%% 
 
def createMetaboliteGraph(m):
    """
    Creates a NetworkX undirected metabolite graph from a Cobrapy metabolic model. In the metabolite graph all 
    metabolites are added as nodes and edges are added between two metabolites if the metabolites occur 
    in the same chemical reaction.

    Parameters
    ----------
    m : Cobrapy metabolic model
    
    Returns
    -------
    Gm : NetworkX metabolite graph 
    dict_size : dictionary with the number of nodes and the number of edges in the graph

    """
    
    #initialize graph 
    Gm = nx.Graph()
    #add all metabolites as nodes to the graph
    for i in range(len(m.metabolites)):
        metabolite = m.metabolites[i]
        Gm.add_node(metabolite.id)
    
    #add all edges to the graph 
    #loop over all metabolites and check in which reactions the metabolite occurs
    for j in range(len(m.metabolites)): 
        metabolite = m.metabolites[j]
        reactions = list(metabolite.reactions) 
        #loop over all reactions in which the metabolite occurs and look which metabolites occur in these reactions
        for r in range(len(reactions)):
            reaction = reactions[r]
            adj_metabolites = list(reaction.metabolites)
            #loop over all metabolites that are connected
            #Add an edge when the metabolite identifiers are not the same. 
            for me in range(len(adj_metabolites)):
                adj_metabolite = adj_metabolites[me]
                
                if metabolite.id != adj_metabolite.id: #do not allow self loops
                    Gm.add_edge(metabolite.id, adj_metabolite.id) 
                    
    nr_nodes = Gm.number_of_nodes()
    nr_edges = Gm.number_of_edges()
    
    dict_size = {"number of nodes": nr_nodes, "number of edges": nr_edges}
                    
    return Gm, dict_size

#%%

def createMetaReacGraph(m):
    """
    Creates a NetworkX undirected metabolite-reaction graph from a Cobrapy metabolic model. In the metabolite graph 
    all metabolites and reactions are added as nodes and edges are added between a metabolite and 
    reaction if the metabolites occurs in the reaction

    Parameters
    ----------
    m : Cobrapy metabolic model

    Returns
    -------
    Gmr : NetworkX undirected metabolite-reaction graph
    dict_size : dictionary with the number of nodes and the number of edges in the graph

    """
    
    #initialize graph 
    Gmr = nx.Graph()
    #add all nodes to the graph
    #add all metabolites as nodes to the graph
    for i in range(len(m.metabolites)):
        metabolite = m.metabolites[i]
        Gmr.add_node(metabolite.id)
    #add all reactions as nodes to the graph
    for j in range(len(m.reactions)):
        reaction =  m.reactions[j]
        Gmr.add_node(reaction.id)
        
    #add all edges
    for k in range(len(m.metabolites)):
        metabolite = m.metabolites[k]
        reactions = list(metabolite.reactions)
        for r in range(len(reactions)):
            reaction = reactions[r]
            Gmr.add_edge(metabolite.id, reaction.id) 

    nr_nodes = Gmr.number_of_nodes()
    nr_edges = Gmr.number_of_edges()
    
    dict_size = {"number of nodes": nr_nodes, "number of edges": nr_edges}
    
    return Gmr, dict_size
              
#%%
def createWeightedGraph(m):
    """
    Creates a weighted and directed metabolite-reaction graph. Edges are only pointing in directions
    corresponding with the reversibility of the reaction.

    Parameters
    ----------
    m : Cobrapy metabolic model

    Returns
    -------
    Gmrw : NetworkX weighted and directed metabolite-reaction graph

    """
    
    s = m.optimize()
    fluxes = s.fluxes
    #initialize graph 
    Gmrw = nx.DiGraph()
    #add all nodes to the graph
    #add all metabolites as nodes to the graph
    for i in range(len(m.metabolites)):
        metabolite = m.metabolites[i]
        Gmrw.add_node(metabolite.id)
    #add all reactions as nodes to the graph
    for j in range(len(m.reactions)):
        reaction =  m.reactions[j]
        Gmrw.add_node(reaction.id)
      
    #add all edges
    for reac in m.reactions:
        metadict = reac.metabolites
        #check if reaction is reversible 
        if reac.lower_bound < 0: #reversible reaction
            #if reaction is reversible add edges in both directions
            for metabolite in metadict:   
                Gmrw.add_edge(reac.id, metabolite.id, weight = abs(fluxes[reac.id]/2))
                Gmrw.add_edge(metabolite.id, reac.id, weight = abs(fluxes[reac.id]/2))
        else:
            #if reaction is irriversible look at the stoichiometric coefficient
            for metabolite in metadict:
                if metadict[metabolite] < 0: 
                    Gmrw.add_edge(metabolite.id, reac.id, weight = abs(fluxes[reac.id]/2))
                else:
                    Gmrw.add_edge(reac.id, metabolite.id, weight = abs(fluxes[reac.id]/2))
                    
    return Gmrw



    