%Script that plots bar graphs with the means senstivity indices per model
%parameters for each of the four state variables looked at (figures in
%section 2 of the report 

%%
%load all data in 
all_SI_G = load("SI_G.mat");
all_SI_I = load("SI_I.mat");
all_SI_NEFA = load("SI_NEFA.mat");
all_SI_TG = load("SI_TG.mat");

%get the sensitivity indices for all parameters
SI_G = all_SI_G.all_SI_tot;
SI_I = all_SI_I.all_SI_tot;
SI_NEFA = all_SI_NEFA.all_SI_tot;
SI_TG = all_SI_TG.all_SI_tot;

%% make plots 
SIs = [SI_G; SI_I; SI_NEFA; SI_TG];
sv_names = {'plasma glucose', 'plasma insulin','plasma NEFA', 'plasma TG'};

for i = 1:4 % for loop over the four state variables in plasma. 
    figure()
    X = categorical({'k_1','k_2','k_3','k_4','k_5','k_6','k_7','k_8','k_9','k_{10}','sigma','K_M','G_b','I_b','EGP_b','f_{spill}','k_{11}','ATL_{max}','K_{ATL}','k_{12}','tau_{LPL}','k_{13}','k_{14}','k_{15}','k_{16}'});
    X = reordercats(X,{'k_1','k_2','k_3','k_4','k_5','k_6','k_7','k_8','k_9','k_{10}','sigma','K_M','G_b','I_b','EGP_b','f_{spill}','k_{11}','ATL_{max}','K_{ATL}','k_{12}','tau_{LPL}','k_{13}','k_{14}','k_{15}','k_{16}'});
    Y = SIs(i,:);
    bar(X,Y)
    ylabel("sensitivity index", 'FontSize',15)
    xlabel('model parameter', 'FontSize',15)
    %title("senstivity indices for model parameters after LSA for state variable plasma", sv_names(i), 'Fontsize', 15)
    set(gca,'FontSize',14);
end
