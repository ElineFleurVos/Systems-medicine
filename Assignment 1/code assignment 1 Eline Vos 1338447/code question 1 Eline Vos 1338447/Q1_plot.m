%This script plots the simulated the Mixed Meal Model for different values
%of the model parameters Km and the sample person characteristic body
%weight (BW). It plots the figures which can be found in section 1 of the
%report. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Specify phenotypic traits necessary for model  simulation
sample_person.glucose = 5;    %fasting glucose (mmol/l)
sample_person.insulin = 18;   %fasting insulin (uIU/ml)
sample_person.TG      = 1.3;  %fasting plasma triglyceride (mmol/l)
sample_person.NEFA    = 0.33; %fasting plasma NEFA(mmol/l)
sample_person.BW      = 84.2; %body weight (kg)

%% specify composition of meal being simulated
sample_person.meal.G  = 75000; %mass of glucose in meal (mg)
sample_person.meal.TG = 60000; %mass of triglyceride in meal (mg)

%% specify value for each model parameter
%glucose + insulin parameters (Rozendaal et al. (2018))
parameters(1) = 0.0105;    %k1 rate constant for glucose stomach emptying (fast)[1/min]
parameters(2) = 0.28;      %k2 rate constant for glucose appearence from gut [1/min]
parameters(3) = 6.07e-3;   %k3 rate constant for suppresstion of hepatic glucose release by change of plasma glucose
parameters(4) = 2.35e-4;   %k4 rate constant for suppression of hepatic glucose release by delayed (remote) insulin
parameters(5) = 0.0424;    %k5 rate constant for delayed insulin depedent uptake of glucose
parameters(6) = 2.2975;    %k6 rate constant for stimulation of insulin production by the change of plasma glucose concentration (beta cell funtion)
parameters(7) = 1.15;      %k7 rate constant for integral of glucose on insulin production (beta cell function)
parameters(8) = 7.27;      %k8 rate constant for the simulation of insulin production by the rate of change in plasma glucose concentration (beta cell function)
parameters(9) = 3.83e-2;   %k9 rate constant for outflow of insulin from plasma to interstitial space
parameters(10) = 2.84e-1;  %k10 rate constant for degredation of insulin in remote compartment
parameters(11) = 1.4;      %sigma shape factor (appearance of meal)
parameters(12) = 13.2;     %Km michaelis-menten coefficient for glucose uptake
parameters(13) = sample_person.glucose(1);%G_b basal plasma glucose [mmol/l]
parameters(14) = sample_person.insulin(1); %I_PL/_b basal plasma glucose [microU/ml]
parameters(15) = 0.043;    %EGP_bbasal hepatic glucose release

%triglyceride + NEFA parameters (new)
parameters(16) = 30;       %f_spill - fractional spillover of LPL derived NEFA
parameters(17) = 0.00045;  %k11 - rate coeficient LPL lipolysis sips/jelic models
parameters(18) = 0.215;    %ATL_max maximum rate of ATL lipolysis in adipose tissues
parameters(19) = 0.0385;   %K_ATL michealis menten coeficient for ATL lipolysis of store TG in adipose tissue
parameters(20) = 0.0713;   %k12 rate constanst for uptake of NEFA into tissues (currently insulin indenpendent)
parameters(21) = 208.88;   %tau_LPL time delay for insulin stimulation of LPL lipolysis
parameters(22) = 0.0088;   %k13 - rate constant for stomach emptying TG(very slow)
parameters(23) = 0.0163;   %k14 - rate constant for rate of TG appearance from gut
parameters(24) = 1e-5;     %k15 coefficient for inhibition of TG secretion from liver by insulin
parameters(25) = 0.0119;   %k16 basal secretion of TG from liver

%% general simulation

%generate plot of model simualtion. 
figure(1)
%specify the time span for the model simulation
time = 0:1:500;
Q1_simulate_meal_model(sample_person,parameters,time);

% %% Simulation for different values of bodyweight
% bws = [50,80,110,130,160]; %define different values for body weights
% parameters1 = parameters;
% sample_person1 = sample_person;
% 
% for i = 1:length(bws)
%     sample_person1.BW = bws(i); 
%     
%     %generate plot of model simualtion. 
%     figure(2)
%     Q1_simulate_meal_model(sample_person1,parameters1,time);
%     hold on 
%     legend('50 kg','80 kg','110 kg','130 kg','160 kg', 'Fontsize', 12);
%     %sgtitle('Simulations for different bodyweights')
% 
%     Lgnd = legend('show');
%     Lgnd.Position(1) = 0.09;
%     Lgnd.Position(2) = 0.2;
%     title(Lgnd,'Body weight')
%         
% end

%% Simulation for different values of model parameter k5
k5s = [0.01,0.02,0.04,0.05,0.07]; %define different values for k5
parameters2 = parameters;

for i = 1:length(k5s)
    parameters2(5) = k5s(i); 
    
    %generate plot of model simualtion. 
    figure(3)
    Q1_simulate_meal_model(sample_person,parameters2,time);
    hold on 
    legend('0.01 min^{-1}','0.03 min^{-1}','0.04 min^{-1}','0.05 min^{-1}','0.07 min^{-1}', 'Fontsize', 12);
    %sgtitle('Simulations for different bodyweights')

    Lgnd = legend('show');
    Lgnd.Position(1) = 0.09;
    Lgnd.Position(2) = 0.2;
    title(Lgnd,'k_5')
        
end

%%  Simulation for different values of model parameter ATL_max
ATLs = [0.05,0.1,0.2,0.3,0.35]; %define different values for ATL_max
parameters3 = parameters;

for i = 1:length(ATLs)
    parameters3(19) = ATLs(i); 
    
    %generate plot of model simualtion. 
    figure(4)
    Q1_simulate_meal_model(sample_person,parameters3,time);
    hold on 
    legend('0.05 min^{-1}','0.1 min^{-1}','0.2 min^{-1}','0.3 min^{-1}','0.5 min^{-1}', 'Fontsize', 12);
    %sgtitle('Simulations for different bodyweights')

    Lgnd = legend('show');
    Lgnd.Position(1) = 0.09;
    Lgnd.Position(2) = 0.2;
    title(Lgnd,'ATL_{max}')
        
end

%% simulation for different values of model parameter Km
Kms = [5,9,13,17,21]; %define different values for Km
parameters4 = parameters;

for i = 1:length(Kms)
    parameters4(12)= Kms(i); 
    
    %generate plot of model simualtion. 
    figure(5)
    Q1_simulate_meal_model(sample_person,parameters4,time);
    hold on 
    legend('5 mmol/L','9 mmol/L','13 mmol/L','17 mmol/L','21 mmol/L', 'Fontsize', 12);
    %sgtitle('Simulations for different bodyweights')

    Lgnd = legend('show');
    Lgnd.Position(1) = 0.09;
    Lgnd.Position(2) = 0.2;
    title(Lgnd,'K_m')
        
end
