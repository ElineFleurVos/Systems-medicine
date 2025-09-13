% This file fits the meal model to the sample data. A multi-start
% optimization approach is implemented. 

%% load the data 
data = load("sample_data.mat");
sample_data = data.sample_data;
 
%% initialize sample person characteristics and initial parameter values
%take the values at timepoint zero as initial values of the state variables
sample_person.glucose = sample_data.glucose(1);
sample_person.insulin = sample_data.insulin(1);   %fasting insulin (uIU/ml)
sample_person.TG      = sample_data.TG(1);  %fasting plasma triglyceride (mmol/l)
sample_person.NEFA    = sample_data.NEFA(1); %fasting plasma NEFA(mmol/l)
%get body weight from sample data 
sample_person.BW      = sample_data.BW;
%get meal information from sample data 
sample_person.meal.G  = 75000; %mass of glucose in meal (mg)
sample_person.meal.TG = 60000; %mass of triglyceride in meal (mg)

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

%% Define subset of model parameters to update during optimization process. 
%SET: put the variables you want to look at here
param = [par0(1), par0(2), par0(3), par0(5), par0(6), par0(7), par0(8), par0(9), par0(10), par0(11), par0(12), par0(13), par0(14), par0(15), par0(18), par0(19), par0(21), par0(22), par0(24), par0(25)]; 
no_param = length(param); %get number of parameters to look at

%% Monte Carlo Multiple Minimization (MCMM) to detect multiple local minima
N = 100; %number of initial parameter sets, must be even number
initial_ensembles_up = zeros(N/2,no_param);
initial_ensembles_down = zeros(N/2,no_param);
scale_up = lhsdesign(N/2,no_param); %latin hybercube sampling
%half of the samples are scaled up and half of samples are scaled down
for i = 1:N/2
    for j = 1:no_param
            initial_ensembles_up(i,j) = param(j) + 0.5*param(j)*scale_up(i,j);
    end
end
scale_down = lhsdesign(N/2,no_param);
for i = 1:N/2
    for j = 1:no_param
            initial_ensembles_down(i,j) = param(j) - 0.5*param(j)*scale_down(i,j);
    end
end
%put initial ensembles scaled up and down together
initial_ensembles = [initial_ensembles_up; initial_ensembles_down];

%% get best parameter sets based on initial errors

time = sample_data.time_G;
initial_errors = zeros(N,1);
for i = 1:N
    %calculate the inital error with meal model error funciton. 
    %If an error occurs, go on with the for loop. 
    try 
        ksi = M3al_Model_ErrorFunc(initial_ensembles(i,:), sample_data, time);
    catch 
        continue
    end
    SSE = sum(ksi.^2); %calculate SSE from differences between true and predicted model outputs
    initial_errors(i) = SSE;
end

% keep best parameter sets
perc = 0.2; %SET: percentage of total amount of parameter sets that you want to keep
nr_keep = perc*N; 
initial_ensembles(:,(length(param)+1)) = initial_errors; %add column with erros after paramaters
[~,idx] = sort(initial_ensembles(:,no_param+1)); % sort the fourth column with initial errors, but NOT zeros in front
ensembles_sorted = initial_ensembles(idx,:); %sort the whole matrix
k = 1;
%now remove the initial values of 0 as these were failed simulations
for i = 1:N
    if ensembles_sorted(i,no_param+1) ~= 0 
        sorted_without_0(k,:) = ensembles_sorted(i,:);
        k = k+1;
    end
end
final_ensembles = sorted_without_0(1:nr_keep,1:no_param); %get the parameter values of best parameter sets

%% fit model to different initial parameter values

%set lowerbounds and upperbounds 
lb = zeros(1,no_param);
ub = [];

%initialize some cells and matrices to store results
all_par_opt = zeros(length(final_ensembles(:,1)), length(param));
final_errors = zeros(1, length(final_ensembles(:,1)));
Mresidual = cell(length(final_ensembles(:,1)), 1);
Mjacobian = cell(length(final_ensembles(:,1)), 1);
parfor i = 1:length(final_ensembles(:,1)) %loop over amount of final ensembles
    %optimize parameters with lsqnonlin using different sets of initial
    %parameter values
    %if error occurs, go on with for loop 
    try 
        result = Fit_M3al_Model(final_ensembles(i,:), sample_data, time, lb, ub);
    catch 
        continue
    end

    %get all the results
    par_opt = result.p_opt;
    resnorm = result.resnorm;
    jacobian = result.jacobian;
    residual = result.residual;

    %put results in cells and matrices for later use.
    final_errors(i) = resnorm;
    all_par_opt(i,:) = par_opt;
    Mresidual{i} = residual;
    Mjacobian{i} = jacobian;
end

%save results to calculate CI and make plots later
all_par_opt(:,no_param+1) = final_errors;
save("opt_params_with_errors2", "all_par_opt")
save("residuals2", "Mresidual")
save("jacobians2", "Mjacobian")


