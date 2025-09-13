# Systems-medicine

This repository contains all code belonging to the course Systems medicine. The course was separated into two assigments.


### Assignment 1: Analysis of mixed meal model

Each question has its own folder which contains the code for this specific question. Zip files of these folders are uploaded separately 
for each of the questions in canvas. Here, for each of the questions it is explained what files are made and what they do. All the 
figures in the report are made with files that have a name of the form QX.plot. All the figures can also be found in the folder named Figures.

QUESTION 1
- Q1_simulate_meal_model: Has some small changes with respect to the file provided with the assignment. Plot colour is now not a parameter of the function anymore.
- Q1_plot: Plots the figures in section 1 of the report.

QUESTION 2
- Q2_LSA: Calculates the mean sensitivity indices per parameter for each state variable. The state variable needs to be specified in line 64 of the file. The series with the mean sensitivity indices are then saved inside the folder. These are already in the folder, so this file does NOT need to be run again. 
- Q2_simulate_meal_model: Has some changes with respect to the file provided with the assignment. It now returns the values of the state variables and does not plot any figures. 
- Q2_plot: Plots the figures in section 2 and Appendix A of the report. 

QUESTION 3(/4)
- Q3_MPSA_with_sp: Function that does the multi-parametric sensitivity analysis. This function also varies the sample person characteristics. The function is used in Q3_plots_MPSA.m.
- Q3_MPSA_without_sp: Function that does the multi-parametric sensitivity analysis. This function does not vary the sample person characteristics. The function is used in Q3_plots_MPSA.m.
- Q3_MPSA_no_sp: This file runs the MPSA analysis for the input you choose. This file is used to run the analysis for large amounts of ensembles, because this took some time to run. 
- Q3_simulate: Here the feature is calculated for a certain ensembled parameter set. This function is used in the MPSA analysis. 
- Q4_choose_parameters: Here calculations are made to see which parameter sensitivities exceed a certain threshold value. This is used to select a subset of parameters.
- Q3_plot: Plots all the figures in section 3 and Appendix B of the report.

QUESTION 5
- Q5_multiple optimization: Here the multi-start optimization approach is implemented and the model is fitted to the sample data. Some results are saved for further use in Q5_plot.m, so this file does NOT need to be run again.
- Q5_meal_model_parameters: Has some changes with respect to the file provided with the assignment, so that the subset of selected parameters is optimized during the optimization process.
- Q5_simulate_meal_model: Has some changes with respect to the file provided with the assignment. The sample data is now also plotted in the figure.
- Q5_plot: Plots the figure and makes the table that can be found in section 5 of the report. 

QUESTION 6
- Q6_simulate_meal_model: Has some changes with respect to the file provided with the assignment. The data points in the new sample datasets are now also plotted as black stars. 
- Q6_bootstrap: Here the new bootstrapped samples are made and the model is fit to these new sample datasets. The results are saved, so this file does not need to be run again, as it takes quite long.
- Q6_plot: Plots the figures and makes the table that can be found in section 6 of the report. 

QUESTION 7
- Q7_mock_datasets: In this file the mock datasets are created by resampling from the fitted model plus adding 10% noise to the data.
- Q7_meal_model_error_func: Has some changes with respect to the file provided with the assignment. Changes are made to the difference between true and predicted model outcomes for state variables NEFA plasma en TG plasma as explained in section 7 of the report. 
- Q7_meal_model_parameters: Has some changes with respect to the file provided with the assignment and the one used in question 5, so that the new subset of selected parameters is optimized during the optimization process.
- Q7_simulate_meal_model: Has some changes with respect to the file provided with the assignment. It plots the simulated data and timepoints for all three sampling schedules. This file also contains calculations of the SSE as found in Table 3 of the report. 
- Q7_run_simulation: Here the model is fitted to the data from the new sample schedule (one sample schedule at a time). A multi-start optimization approach is used. 
- Q7_plot: Plots the figure that can be found in section 7 of the report. 

### Assignment 2 Bioinformatics

This read me contains information on the files and code which are found in the zip file. The zip file contains of four folders:
'Code part A', 'Code part B', 'Bacteria' and 'Figures'. The latter contains all the figures that can be found in the report. 'Bacteria'
contains the bacteria sbml files used in the study. The first two folders contain the following code files:

CODE PART A:
- ga_algorithm: Self-designed global alignment algorithm.
- load_sequence: Loads the reference sequence (or a any sequence in general) from the NCBI library. 
- check_chunks: This file checks the difference between calculating the similarity score as a sum of two chunked sequences
                compared to doing the semi global alignment in one go. Also,the similarity score and total 
                amount of matches is calculated when using the maximal chunk size of 5000. 
- analysis_proteins: Does semi-global alignment for the different protein coding sequences and calculates the number of matches and 
                     percentage of nucleotides with a match. 
- translation_function: Contains a function that translate a RNA sequene to an amino acids sequence.
- analysis_spike: Analyzes the spike protein by finding the exact match of the RNA sequences and translating these RNA sequences to 
                  amino acid sequences.
- analysis_membrane_protein: Analyses the membrane protein semi-global alignment
- analysis_ORF_7a_8: Checks if the ORF7a and ORF8 mutations of Delta 21A are actually present on the reference sequence NC_045512. 

- score_matrix.xlsx: Excel file that contains the score matrix as found in the report. Colors are added here manually. 

CODE PART B:
- create_graphs: Contains functions that make an undirected metabolite graph, undirected metabolite-reaction graph and directed, weighted
                 metabolite-reaction graph (from a Cobrapy model).
- plot_graphs: Plots the graphs found in part B of the report. 
- static_analysis: Plots the rank-degree plots, calculates hubs and plots the 10 metabolites (hubs) with highest degree.
- shortest_path: Contains one function that calculates the shortest paths from metabolites in the medium to the biomass and two functions
                 necessary to plot bar plots with the most frequent lengths of the shortest paths.  
- change_media: Contains all functions that are associated with the FBA in DMEM and an empty medium. 
- FBA: Runs and plots the results of the FBA in the medium defined by the SBML file.  
- FBA_DMEM_empty: Runs and plots the results of the FBA in DMEM and the empty medium.
- FBA_symbiosis: Contains all code made so far to look at interaction between bacteria. 
