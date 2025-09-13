% this file runs the optimization procedure for ONE sample strategy including a multistart
% optimization approach. If you want results for all three sample strategy,
% you need to run this file three times. 

%% load data and initial parameters
par0 = load("parameters.mat").par0;
mock_data = load("mock_test.mat").mock_data;

sample_data0 = mock_data(1);
sample_data1 = mock_data(2);
sample_data2 = mock_data(3);
sample_data3 = mock_data(4);

%CHOOSE SAMPLE SET YOU WANT TO LOOK AT
sample_data = sample_data3;
name = "result_sstest";

%% initialize sample person characteristics 
%take the values at timepoint zero as initial values of the state variables
sample_person.glucose = sample_data.glucose(1);
sample_person.insulin = sample_data.insulin(1);   %fasting insulin (uIU/ml)
sample_person.TG      = sample_data.TG(1);  %fasting plasma triglyceride (mmol/l)
sample_person.NEFA    = sample_data.NEFA(1); %fasting plasma NEFA(mmol/l)
%get body weight from sample data 
sample_person.BW      = sample_data.BW;
%get meal information from sample data 
sample_person.meal.G  = sample_data.meal.G; %mass of glucose in meal (mg);

%% Define subset of model parameters to update during optimization process. 

%define the trainable parameters
param = [par0(1), par0(5), par0(11), par0(13), par0(14), par0(18), par0(19), par0(21), par0(22)]; 
no_param = length(param); %number of trainable parameters

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
    %if an error occurs, go on with the for loop. 
    try 
        ksi = Q7_meal_model_error_func(initial_ensembles(i,:), sample_data, time);
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

final_errors = zeros(1, length(final_ensembles(:,1)));
all_par_opt = zeros(length(final_ensembles(:,1)), length(param));
parfor i = 1:length(final_ensembles(:,1)) %loop over amount of final ensembles
    %optimize parameters with lsqnonlin using different sets of initial
    %parameter values
    %if error occurs, go on with for loop 
    try 
        result = Fit_M3al_Model(final_ensembles(i,:), sample_data, time, lb, ub);
    catch 
        continue
    end
    %get optimum parameters and ksi
    par_opt = result.p_opt;
    resnorm = result.resnorm;
    
    final_errors(i) = resnorm;
    all_par_opt(i,:) = par_opt;
end
all_par_opt(:,no_param+1) = final_errors; %add SSEs as new most right column

%save data 
save(name, "all_par_opt");


