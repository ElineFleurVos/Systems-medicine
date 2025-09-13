% this code plots the figure and makes the table found in section 5 of the report. 

%% load results
opt_params_errors0 = load("opt_params_with_errors.mat").all_par_opt;
residuals = load("residuals.mat").Mresidual;
jacobians = load("jacobians.mat").Mjacobian;

% remove failed optimizations
k = 1;
for i = 1:length(opt_params_errors0(:,1))
    if opt_params_errors0(i,1) > 0
        opt_params_errors(k,:) = opt_params_errors0(i,:);
        k = k+1;
    end
end

%% find optimum parameters and index with lowest final error

[~,idx] = sort(opt_params_errors(:,end)); % sort the last column with final errors
errors_sorted = opt_params_errors(idx,:); %sort the whole matrix

%% calculate 95% confidence interval for best paramameters
no_param = length(errors_sorted)-1; %number of parameters

best_idx = idx(1);
best_par_opt = errors_sorted(1,1:no_param);
best_residual = residuals{best_idx};
best_jacobian = jacobians{best_idx};
best_error = errors_sorted(1,end); %get the best error out of the first row and final column

N = 24;
best_covarpar = best_error*inv(best_jacobian'*best_jacobian)/N;
CI = nlparci(best_par_opt,best_residual,'covar', best_covarpar);

parameter = ["k_1";"k_2";"k_3";"k_5";"k_6";"k_7";"k_8";"k_9";"k_{10}";"sigma";"K_M";"G_b";"I_b";"EGP_b";"ATL_{max}";"K_{ATL}";"tau_{LPL}";"k_{13}";"k_{15}";"k_{16}"];
value = transpose(best_par_opt);
low_bound = CI(:,1);
up_bound = CI(:,2);
T = table(value,low_bound,up_bound);
T.Properties.RowNames = parameter;
table2latex(T, 'CI_results');

save("best_par", "best_par_opt")

%% plotting 

data = load("sample_data.mat");
sample_data = data.sample_data;

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

save("parameters", "par0")

all_opt_params = opt_params_errors(:,1:no_param); %remove the columns with the errors

%specify plotting colour
plot_colour = [0, 0.4470, 0.7410]; 
%generate plot of model simulation. 
figure(1)
%specify the time span for the model simulation
time = 0:1:500;

nr_lines = 5; %SET: choose the number of simulated lines in the plot. 
num_params = [1,2,3,5,6,7,8,9,10,11,12,13,14,15,18,19,21,22,24,25];
for s = 1:nr_lines
    opt_par = all_opt_params(s,:);
    par_new = par0;
    for i = 1:no_param
        num = num_params(i);
        par_new(num) = opt_par(i);
    end
    Q5_simulate_meal_model(sample_person,par_new,time,plot_colour);
    hold on
end


