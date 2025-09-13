%Script to perform local parameter sensitivity analysis for the Mixed Meal Model. 
% 
% This script calculates the sensitivity indices for all 25 model
% parameters for different state variables. These sensitivity indices are
% then saved and plots are made in matlab file ex2_LSA.plot.m. 

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
parameters(1) = 0.0105;  %k1 rate constant for glucose stomach emptying (fast)[1/min]
parameters(2) = 0.28;      %k2 rate constant for glucose appearence from gut [1/min]
parameters(3) = 6.07e-3;   %k3 rate constant for suppresstion of hepatic glucose release by change of plasma glucose
parameters(4) = 2.35e-4;   %k4 rate constant for suppression of hepatic glucose release by delayed (remote) insulin
parameters(5) = 0.0424;  %k5 rate constant for delayed insulin depedent uptake of glucose
parameters(6) = 2.2975;  %k6 rate constant for stimulation of insulin production by the change of plasma glucose concentration (beta cell funtion)
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
parameters(17) = 0.00045; %k11 - rate coeficient LPL lipolysis sips/jelic models
parameters(18) = 0.215;    %ATL_max maximum rate of ATL lipolysis in adipose tissues
parameters(19) = 0.0385; %K_ATL michealis menten coeficient for ATL lipolysis of store TG in adipose tissue
parameters(20) = 0.0713; %k12 rate constanst for uptake of NEFA into tissues (currently insulin indenpendent)
parameters(21) = 208.88; %tau_LPL time delay for insulin stimulation of LPL lipolysis
parameters(22) = 0.0088;   %k13 - rate constant for stomach emptying TG(very slow)
parameters(23) = 0.0163; %k14 - rate constant for rate of TG appearance from gut
parameters(24) = 1e-5;     %k15 coefficient for inhibition of TG secretion from liver by insulin
parameters(25) = 0.0119; %k16 basal secretion of TG from liver

%specify range through which each parameter will be varied
range = 0.5; %50% 
%specify number of parameter values in range to be simulated. 
n=6; 
%specify the time span for the model simulation
t_span = 0:1:500;
%GIVE THE STATE VARIABLE YOU WANT TO LOOK AT AS INPUT HERE 
sv = 13; %% choose which state variable to look at, G_PL = 2, I_PL = 4, NEFA_PL = 9, TG_PL = 13

%preallocate mean sensitivity indices
all_SI_tot = zeros(1,length(parameters));

for p = 1:length(parameters) %for loop over all 25 model parameters

    [T, X] = Q2_simulate_meal_model(sample_person,parameters,t_span);
    G_PL = X(:,sv);
    
    %p 
    p_reference = parameters(p);             %define reference value
    p_lb=parameters(p)-range*parameters(p);  %define lower bound parameter
    p_ub=parameters(p)+range*parameters(p);  %define upper bound parameter
    step=(p_ub-p_lb)/(n-1);
    p_values=p_lb:step:p_ub;                 %define all n parameter values to be tested

    parameter_values=repmat(parameters,n,1);
    parameter_values(:,p)=p_values';

    SI_p = zeros(1,length(p_values)); %pre-allocate all sensitivity indices
    for i = 1:length(p_values) %for loop over the different new parameter values 
        p_new = p_values(i);
        delta_p = abs(p_new-p_reference);
        % set try catch if error is found
        try [T_new, X_new] = Q2_simulate_meal_model(sample_person,parameter_values(i,:),t_span);
        catch e
            continue
        end
        G_PL_new = X_new(:,sv);
    
        SIt = zeros(1,length(G_PL)); %pre-allocate SIs for all time points
        for j = 1:length(G_PL) %for loop over all time points 
                delta_M = abs(G_PL_new(j) - G_PL(j));
                M0 = G_PL(j);
                SI_j = (delta_M/delta_p)*(p_reference/M0); %calculate sensitivity index per time point
                SIt(j) = SI_j;
        end
        SI_p(i) = sum(SIt); %calculate the sum of the sensitivity indices over all time points. 
    end

    SI_tot = sum(SI_p)/length(p_values); %calculate the mean SI over all different parameter values
    all_SI_tot(p) = SI_tot;
end

save('SI_test.mat',"all_SI_tot")
    