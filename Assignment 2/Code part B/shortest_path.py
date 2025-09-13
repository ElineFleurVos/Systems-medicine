# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 22:33:32 2022

@author: 20182460
"""
import cobra
import matplotlib.pyplot as plt
import networkx as nx
from create_graphs import createWeightedGraph
import operator

#%%
def shortestPath(m):
    """
    Calculates the length of the shortest path from a metabolite in the medium to the biomass using the
    Dijkstra algorithm. The length of the path is defined as the sum of the absolute fluxes on the path. 
    The maximum and minimum shortest path length as well as the 3 path lenths with the highest frequency are found.

    Parameters
    ----------
    m : cobrapy metabolic model

    Returns
    -------
    min_value (float): Minimum shortest path found
    max_value (float) : Maximum shortest path found
    list_lengths_str (list): contains the 3 path lengths with the highest frequency
    list_freq (list ): contains the frequences of the 3 path lengths with the highest freqency

    """
    
    #create weighted graph
    Gmrw = createWeightedGraph(m)
    #get medium
    medium = m.medium
    path_lengths = {}
    list_lengths = []
    medium_names = list(medium.keys())
    for i in range(len(medium_names)): #loop over metabolites in the medium
        metabolite = medium_names[i]
        try:
            #calculate shortest path length using Dijkstra algorithm
            length = nx.dijkstra_path_length(Gmrw, metabolite, target='biomass_c', weight='weight')
        except: 
            continue
            
        path_lengths[metabolite] = length #add path length to dictionary
        list_lengths.append(length)
        
    #get the metabolite name and shortest path lenght of the minimum and maximum shortest path lengths.
    min_key = min(path_lengths, key=path_lengths.get)
    min_value = path_lengths[min_key]
    max_key = max(path_lengths, key=path_lengths.get)
    max_value = path_lengths[max_key]
    
    #calculate for each unique shortes path length the frequency
    frequency_dict = {}
    unique_values = list(set(path_lengths.values())) #find all unique path lengths
    for i in range(len(unique_values)): #loop over unique path lengts
        frequency = 0
        for j in list(path_lengths.keys()): #loop over all metabolies
            if path_lengths[j] == unique_values[i]:
                frequency += 1
        frequency_dict[unique_values[i]] = frequency
        
    #shrink dictionary to contain only 3 path lengts with highest frequency
    high_freq_dict = dict(sorted(frequency_dict.items(), key=operator.itemgetter(1), reverse=True)[:3])
        
    #make lists of the key of the dictionary and the values of the dictionary for plotting later.
    list_lengths = list(high_freq_dict.keys())
    list_lengths_str= []
    list_freq = []
    for i in range(len(list_lengths)):
        rounded = round(list_lengths[i], 1)
        list_lengths_str.append(str(rounded))
        list_freq.append(high_freq_dict[list_lengths[i]])

    return min_value, max_value, list_lengths_str, list_freq


def plotMostFrequent(m1, m2, m3):
    """
    plots bar plots with for each of the bacteria the 3 path lengths with the highest frequency.

    Parameters
    ----------
    m1 : Cobrapy metabolic model
    m2 : Cobrapy metabolic model
    m3 : Cobrapy metabolic model

    Returns
    -------
    ax : matplotlib.axes.Axes object of the figure

    """
    min_value1, max_value1, list_lengths1, list_freq1 = shortestPath(m1)
    min_value2, max_value2, list_lengths2, list_freq2 = shortestPath(m2)
    min_value3, max_value3, list_lengths3, list_freq3 = shortestPath(m3)
    
    total_list_lengths = list_lengths1 + list_lengths2 + list_lengths3
    total_list_freq = list_freq1 + list_freq2 + list_freq3
    
    colors = {'Aero1':'royalblue', 'Aero2':'red', 'Clos':'orange'}         
    labels = list(colors.keys())
    handles = [plt.Rectangle((0,0),1,1, color=colors[label]) for label in labels]
    bar_colors = ['tab:blue','tab:blue','tab:blue','tab:red','tab:red','tab:red','tab:orange','tab:orange','tab:orange']
    
    fig, ax = plt.subplots()
    ax.bar(total_list_lengths, total_list_freq, color=bar_colors)
    ax.set_xlabel('length shortest path', fontsize =  12)
    ax.set_ylabel('frequency',fontsize =  12)
    ax.legend(handles, labels, title = 'Bacterium')
    
    return ax


def add_value_labels(ax, spacing=5):
    """Add labels to the end of each bar in a bar chart.

    Parameters:
    ----------
    ax (matplotlib.axes.Axes) : The matplotlib object containing the axes of the plot to annotate.
    spacing (int): The distance between the labels and the bars.
        
    Returns: None
    """

    # For each bar: Place a label
    for rect in ax.patches:
        # Get X and Y placement of label from rect.
        y_value = rect.get_height()
        x_value = rect.get_x() + rect.get_width() / 2

        # Number of points between bar and label. Change to your liking.
        space = spacing
        # Vertical alignment for positive values
        va = 'bottom'

        # If value of bar is negative: Place label below bar
        if y_value < 0:
            # Invert space to place label below
            space *= -1
            # Vertically align label at top
            va = 'top'

        # Use Y value as label and format number with one decimal place
        label = "{:.1f}".format(y_value)

        # Create annotation
        ax.annotate(
            label,                      # Use `label` as label
            (x_value, y_value),         # Place label at end of the bar
            xytext=(0, space),          # Vertically shift label by `space`
            textcoords="offset points", # Interpret `xytext` as offset in points
            ha='center',                # Horizontally center label
            va=va)                      # Vertically align label differently for
                                        # positive and negative values.