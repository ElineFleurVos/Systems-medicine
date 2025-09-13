% This file plots all the figures which are found in section 6 of the report.

%% Specify phenotypic traits necessary for model  simulation
sample_person.glucose = 5;    %fasting glucose (mmol/l)
sample_person.insulin = 18;   %fasting insulin (uIU/ml)
sample_person.TG      = 1.3;  %fasting plasma triglyceride (mmol/l)
sample_person.NEFA    = 0.33; %fasting plasma NEFA(mmol/l)
sample_person.BW      = 84.2; %body weight (kg)
% specify composition of meal being simulated
sample_person.meal.G  = 75000; %mass of glucose in meal (mg)
sample_person.meal.TG = 60000; %mass of triglyceride in meal (mg)

%load the result in 
par0_all = load('parameters.mat').par0; %load the initial parameter values
bootstrap_data = load("bootstrap_data.mat").all_data;
bootstrap_par0 = load("bootstrap_par").all_par_opt; 

num_params = [1,2,3,5,6,7,8,9,10,11,12,13,14,15,18,19,21,22,24,25]; %trainable parameters
no_params = length(num_params); %number of parameters
par0 = par0_all(num_params); %get initial parameter values for initial parameters

%remove rows with only zeros, which are failed simulations
k = 1;
for i = 1:length(bootstrap_par0)
    if bootstrap_par0(i,:) ~= 0 
        bootstrap_par(k,:) = bootstrap_par0(i,:);
        k = k+1;
    end
end

%% PLotting simulations

%generate plot of model simulation
figure(1)
%specify the time span for the model simulation
time = 0:1:500;
%specify time points
TS = bootstrap_data(1).time_G;

%speficiy colors. One for the true model and one for the simulated models
%from new parameter sets.
plot_col = [0.8500 0.3250 0.0980];
plot_col_true = [0,0,0];

%reverse the number of samples, so that the first row in bootstrap_par is
%the last in the for loop. Then the black line will be in front of all the
%orange lines. Bit cumbersome, but it works.
samples = 1:length(bootstrap_par(:,1));
sample_reverse = flip(samples);

for s = 1:length(bootstrap_par(:,1)) %loop over the number of sample datasets
    ss = sample_reverse(s);
    
    %specify plot colour. The first row in bootstrap_par are the parameters
    %found when fitted to the true data. The line simulated from these
    %parameters needs to be black.
    if ss == 1
        col = plot_col_true;
    else 
        col = plot_col;
    end
    
    %Add the trainable parameters into the total set of parameters.
    par_opt = bootstrap_par(s,:);
    all_par = par0_all;
    for i = 1:no_params
        num = num_params(i);
        all_par(num) = par_opt(i);
    end
    bootstrap = bootstrap_data(s);
    
    G = bootstrap.glucose';
    I = bootstrap.insulin';
    TG = bootstrap.TG';
    NEFA = bootstrap.NEFA';
    bootstrap_sol = [G,I,TG,NEFA];

    Q6_simulate_meal_model(sample_person, all_par, bootstrap_sol, TS, time, col)
    hold on
end

%% Plotting boxplots

figure(2)
CI = zeros(20,2);

parameter = ['k_1';"k_2";"k_3";"k_5";"k_6";"k_7";"k_8";"k_9";"k_{10}";"sigma";"K_M";"G_b";"I_b";"EGP_b";"ATL_{max}";"K_{ATL}";"tau_{LPL}";"k_{13}";"k_{15}";"k_{16}"];
for i = 1:no_params
    par = bootstrap_par(:,i);

    CI_low = prctile(par, 2.5); %calculate lower bound
    CI(i,1) = CI_low;
    CI_up = prctile(par, 97.5); %calculate upper bound
    CI(i,2) = CI_up;
    L = CI_up - CI_low;
    perL = 100*(L/par0(i)); %calculate %L (see report for explanation)
    CI(i,3) = perL;

    subplot(2,10,i) 
    boxplot(par)
    title(parameter(i))
end

lower_bound = CI(:,1);
upper_bound = CI(:,2);
perL = CI(:,3);
parameter_names = ["k_1";"k_2";"k_3";"k_5";"k_6";"k_7";"k_8";"k_9";"k_{10}";"sigma";"K_M";"G_b";"I_b";"EGP_b";"ATL_{max}";"K_{ATL}";"tau_{LPL}";"k_{13}";"k_{15}";"k_{16}"];
T = table(lower_bound, upper_bound, perL); 
T.Properties.RowNames = parameter_names;
table2latex(T, 'CI_results_bootstrap');