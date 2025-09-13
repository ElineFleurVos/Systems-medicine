function [Np, K_S, mean_dummies] = Q3_MPSA_with_sp(nr_dummies, nr_MPSA, sign)
%Does multi-Parametric Sensitivity Analysis with Try-Catch. All 25 model
%parameters as well as the sample person characteristics are varied in the ensembling of 
%parameter sets.

%input
%nr_dummies: the number of dummies used in the multi-parametric sensitivity analysis (integer)
%nr_MPSA: the number of parameters ensembles (integer)
%sign: 1 if you want to add the scaling to the initial parameter values
%and 0 if you want to substract the scaling of the initial parameter values

%output
%Np: total number of parameters (integer)
%K_S: Kolmogorov-Smirnov statistic for all parameters and dummies (double series)
%mean_dummies: mean Kolmogorov-Smirnov statistic for dummies (float)


%% Set reference model parameters used to simulate model and calculate feature  
% Specify phenotypic traits necessary for model  simulation
sample_person.glucose = 5;    %fasting glucose (mmol/l)
sample_person.insulin = 18;   %fasting insulin (uIU/ml)
sample_person.TG      = 1.3;  %fasting plasma triglyceride (mmol/l)
sample_person.NEFA    = 0.33; %fasting plasma NEFA(mmol/l)
sample_person.BW      = 84.2; %body weight (kg)

% specify composition of meal being simulated
sample_person.meal.G  = 75000; %mass of glucose in meal (mg)
sample_person.meal.TG = 60000; %mass of triglyceride in meal (mg)

% specify value for each model parameter
% glucose + insulin parameters (Rozendaal et al. (2018))
par0(1) = 0.0105;    %k1 rate constant for glucose stomach emptying (fast)[1/min]
par0(2) = 0.28;      %k2 rate constant for glucose appearence from gut [1/min]
par0(3) = 6.07e-3;   %k3 rate constant for suppresstion of hepatic glucose release by change of plasma glucose
par0(4) = 2.35e-4;   %k4 rate constant for suppression of hepatic glucose release by delayed (remote) insulin
par0(5) = 0.0424;    %k5 rate constant for delayed insulin depedent uptake of glucose
par0(6) = 2.2975;    %k6 rate constant for stimulation of insulin production by the change of plasma glucose concentration (beta cell funtion)
par0(7) = 1.15;      %k7 rate constant for integral of glucose on insulin production (beta cell function)
par0(8) = 7.27;      %k8 rate constant for the simulation of insulin production by the rate of change in plasma glucose concentration (beta cell function)
par0(9) = 3.83e-2;   %k9 rate constant for outflow of insulin from plasma to interstitial space
par0(10) = 2.84e-1;  %k10 rate constant for degredation of insulin in remote compartment
par0(11) = 1.4;      %sigma shape factor (appearance of meal)
par0(12) = 13.2;     %Km michaelis-menten coefficient for glucose uptake
par0(13) = sample_person.glucose(1);%G_b basal plasma glucose [mmol/l]
par0(14) = sample_person.insulin(1); %I_PL/_b basal plasma glucose [microU/ml]
par0(15) = 0.043;    %EGP_bbasal hepatic glucose release

%triglyceride + NEFA parameters (new)
par0(16) = 30;       %f_spill - fractional spillover of LPL derived NEFA
par0(17) = 0.00045;  %k11 - rate coeficient LPL lipolysis sips/jelic models
par0(18) = 0.215;    %ATL_max maximum rate of ATL lipolysis in adipose tissues
par0(19) = 0.0385;   %K_ATL michealis menten coeficient for ATL lipolysis of store TG in adipose tissue
par0(20) = 0.0713;   %k12 rate constanst for uptake of NEFA into tissues (currently insulin indenpendent)
par0(21) = 208.88;   %tau_LPL time delay for insulin stimulation of LPL lipolysis
par0(22) = 0.0088;   %k13 - rate constant for stomach emptying TG(very slow)
par0(23) = 0.0163;   %k14 - rate constant for rate of TG appearance from gut
par0(24) = 1e-5;     %k15 coefficient for inhibition of TG secretion from liver by insulin
par0(25) = 0.0119;   %k16 basal secretion of TG from liver

% Simulate default model for nominal parameter values par0 and calculate feature from selected output(s)
feature0 = Q3_simulate(sample_person, par0);

%add the sample person characteristics to the parameter list, so these can
%be varied too.
par0(26) = sample_person.glucose;
par0(27) = sample_person.insulin;
par0(28) = sample_person.TG;
par0(29) = sample_person.NEFA;
par0(30) = sample_person.BW;
par0(31) = sample_person.meal.G;
par0(32) = sample_person.meal.TG;

Np = length(par0);  %no. of parameters included in analysis
Nd = nr_dummies;    %no. of dummies
n_MPSA = nr_MPSA;   %no. of parameter ensembles / no. of Monte Carlo simulations
tot_params = Np+Nd;

% Latin Hypercube Sampling of parameter space
scale = lhsdesign(n_MPSA, tot_params);   % random uniform distributed parameter sets (values [0,1])
                                         % (scale is an n_MPSA x Np+Nd matrix)

for i = length(par0)+1:length(par0)+1+Nd % set initial dummy parameter values to 1
     par0(i) = 1;
end

% Use the scale matrix to generate MC parameter sets relative to reference values par0:
par_MC = zeros(n_MPSA, tot_params);
for i = 1:n_MPSA
    for j = 1:tot_params
        if sign == 1
            par_MC(i,j) = par0(j) + scale(i,j)*par0(j); %scale with a maximum of 50%
        else 
            par_MC(i,j) = par0(j) - scale(i,j)*par0(j);
        end
    end
end

%% Monte Carlo simulations
for j = 1 : n_MPSA
    j
    partemp=par_MC(j,1:25); %select the model parameters from set j (NOT sample person characteristics)

    %make a struct of the new sample person characteristics so that
    %simulation still works
    sample_person_new.glucose = par_MC(j,26);
    sample_person_new.insulin = par_MC(j,27);
    sample_person_new.TG      = par_MC(j,28);
    sample_person_new.NEFA    = par_MC(j,29);
    sample_person_new.BW      = par_MC(j,30);

    % specify composition of meal being simulated
    sample_person_new.meal.G  = par_MC(j,31);
    sample_person_new.meal.TG = par_MC(j,32);

    % Simulate model for partemp and calculate feature of interest from
    % selected model outputs
    try %if simulation works, the sensitivity criterion for parameter set j gets a value 
        feature = Q3_simulate(sample_person_new, partemp);
        % Calculate sensitivity criterion
        % Sum of differences between the perturbed and reference output
        V(j) = sum(abs(feature0-feature));
    catch %if simulation fails, the sensitivity criterion for parameter set j is marked as NaN
        V(j)= NaN;
    end
end

%% Classification of acceptable vs unacceptable simulations
flag  = zeros(n_MPSA, 1);
Sa    = zeros(n_MPSA, Np);
Su    = zeros(n_MPSA, Np);
value = zeros(n_MPSA, Np);

threshold = nanmean(V);    %mean as threshold
acc       = find(V <= threshold);
unacc     = find(V  > threshold);
flag(acc) = 1;

%% Cumulative distributions (for model parameters and dummies)
for i = 1 : Np+Nd
    temp = [par_MC(:,i), flag];     %associate 1 to acceptable cases and 0 to unacceptable parameter values
    temp = sortrows(temp,1);        %sorts temp based on column 1
    
    value(:,i) = temp(:,1);
    Sa(:,i) = cumsum(temp(:,end));
    Su(:,i) = cumsum(-1*temp(:,end)+ones(n_MPSA,1));
    Sa(:,i) = Sa(:,i)/max(Sa(:,i));
    Su(:,i) = Su(:,i)/max(Su(:,i)); 
end

%% Kolmogorov-Smirnov statistic
K_S = max(abs(Sa-Su));
mean_dummies = mean(K_S(Np+1:end)); %calculate mean of dummies
