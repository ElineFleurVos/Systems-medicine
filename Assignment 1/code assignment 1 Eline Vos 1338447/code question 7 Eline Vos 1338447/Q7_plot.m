% this file plots the figure found in section 7 of the report

%% load data and initial parameters
parameters = load("parameters.mat").par0;
result_ss0 = load("result_ss0.mat").all_par_opt;
result_ss1 = load("result_ss1.mat").all_par_opt;
result_ss2 = load("result_ss2.mat").all_par_opt;
result_ss3 = load("result_ss3.mat").all_par_opt;

%% Specify phenotypic traits and model parameters
sample_person.glucose = 5;    %fasting glucose (mmol/l)
sample_person.insulin = 18;   %fasting insulin (uIU/ml)
sample_person.TG      = 1.3;  %fasting plasma triglyceride (mmol/l)
sample_person.NEFA    = 0.33; %fasting plasma NEFA(mmol/l)
sample_person.BW      = 84.2; %body weight (kg)

% specify composition of meal being simulated
sample_person.meal.G  = 75000; %mass of glucose in meal (mg)
sample_person.meal.TG = 60000; %mass of triglyceride in meal (mg)

%% Find the best final parameters and highest SSE for each sample schedule
num_params = [1,5,11,13,14,18,19,21,22];
no_param = length(num_params);

[~,idx_ss0] = sort(result_ss0(:,end)); % sort the last column with final errors
errors_sorted_ss0 = result_ss0(idx_ss0,:); %sort the whole matrix
opt_SSE_ss0 = errors_sorted_ss0(1,end);
opt_par_ss0 = errors_sorted_ss0(7,1:no_param);

[~,idx_ss1] = sort(result_ss1(:,end)); % sort the last column with final errors
errors_sorted_ss1 = result_ss1(idx_ss1,:); %sort the whole matrix
opt_SSE_ss1 = errors_sorted_ss1(1,end);
opt_par_ss1 = errors_sorted_ss1(1,1:no_param);

[~,idx_ss2] = sort(result_ss2(:,end)); % sort the last column with final errors
errors_sorted_ss2 = result_ss2(idx_ss2,:); %sort the whole matrix
opt_SSE_ss2 = errors_sorted_ss2(1,end);
opt_par_ss2 = errors_sorted_ss2(1,1:no_param);

[~,idx_ss3] = sort(result_ss3(:,end)); % sort the last column with final errors
errors_sorted_ss3 = result_ss3(idx_ss3,:); %sort the whole matrix
opt_SSE_ss3 = errors_sorted_ss3(2,end); %first row only contained zeros because of failed simulation
opt_par_ss3 = errors_sorted_ss3(2,1:no_param);

%% Plotting

all_opt_par = [opt_par_ss0; opt_par_ss1; opt_par_ss2; opt_par_ss3];

for s = 1:4 %loop over the sampling schedules
    opt_par = all_opt_par(s,:);
    par = parameters;
    %add new optimal parameter to complete parameter set
    for i = 1:no_param
        num = num_params(i);
        par(num) = opt_par(i);
    end
    all_par(s,:) = par;
end

figure(1)
time = 0:1:500; %specify time span
errors = Q7_simulate_meal_model(sample_person, all_par, time);
