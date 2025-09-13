# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 10:37:35 2022

@author: 20182460
"""

import cobra
from cobra.medium import minimal_medium

#%%

#load in the data 
#SET CORRECT PATH NAME
path_name = 'C:/Users/20182460/Desktop/Quartile 1/Systems medicine/Assignment 2/Assignment 2 8CM00 Eline Vos/'
m_aero1 = cobra.io.read_sbml_model(path_name + '/Bacteria/Aeromonas_jandaei_Riv2.xml')
m_aero2 = cobra.io.read_sbml_model(path_name + '/Bacteria/Aeromonas_hydrophila_subsp_hydrophila_ATCC_7966.xml')
m_clos = cobra.io.read_sbml_model(path_name + '/Bacteria/Clostridium_perfringens_ATCC_13124.xml')

#%%
#look at individual minimal medium
max_growth = m_aero1.slim_optimize()
min_medium_aero1 = minimal_medium(m_aero1, max_growth)

max_growth = m_aero2.slim_optimize()
min_medium_aero2 = minimal_medium(m_aero2, max_growth)

medium_clos = m_clos.medium
max_growth = m_clos.slim_optimize()
min_medium_clos = minimal_medium(m_clos, max_growth)

#%% look at interaction between m_aero1 and m_aero2

def combineBacteria(m1, m2):
    """
    Add the reactions and metabolites of model m2 to model m1. The names of the metabolites are changed
    so that the bacteria can still be kept apart. This code also 

    Parameters
    ----------
    m1 : Cobrapy metabolic model where another model is added to 
    m2 : Cobrapy metabolic model that is added to the first model
    
    Returns
    -------
    m12 : New Cobrapy metabolic model with both bacteria in a medium

    """
    m12 = m1.copy()
    #add metabolites from Aero2 to m_aero1_aero2
    reactions2 = m2.reactions
    for i in range(len(reactions2)):
        r = reactions2[i]
        r.id = r.id + '_b2' #change the reaction name
        dict_metabolites = r.metabolites
        metabolites = list(dict_metabolites.keys())
        
        #loop over all metabolites occuring in the reaction
        for m in range(len(metabolites)):
            metabolite = metabolites[m]
            name = metabolite.id 
            
            #this part checks if the name of the metabolit is already changes. When this is not the case, 
            #the name is changed.
            if  name.endswith('_b2') == False:
                #metabolite.id = metabolite.id + '_b2'
                coefficient = dict_metabolites[metabolite]
                new_metabolite = metabolite
                new_metabolite.id = metabolite.id + '_b2'
                dict_metabolites[new_metabolite] = coefficient   
                
    #append all the new reaction with new names
    new_reactions = list(m2.reactions)
    for reac in new_reactions:
        m12.add_reaction(reac) 

    return m12

#%%
m_aero1_aero2 = combineBacteria(m_aero1, m_aero2) 

#%% make medium

#get the biomasses of both bacteria
biomass_aero1 = m_aero1_aero2.reactions.get_by_id("biomass423")
biomass_aero2 = m_aero1_aero2.reactions.get_by_id("biomass444_b2")
#set the ojective funtion as the sum of the biomass reaction of both bacteria
m_aero1_aero2.objective = {biomass_aero1:1, biomass_aero2:1}

#get the minium medium in which both bacteria can growth with same sum of fluxes of the 
#biomass reactions as before
max_growth = m_aero1_aero2.slim_optimize()
min_medium = minimal_medium(m_aero1_aero2, max_growth)
m_aero1_aero2.medium = min_medium #set new minimum medium
#check if objective value the same as before change of medium
s_min = m_aero1_aero2.optimize() #optimize
new_score = s_min.objective_value
#new score is indeed still the same

#%% add interactions

summary = m_aero1_aero2.summary()
secretion = summary.secretion_flux
uptake = summary.uptake_flux



            
