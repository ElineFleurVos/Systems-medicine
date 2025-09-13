% This file does the bootstrap sampling en fits the model to the new
% sampled datasets. 

%% loading in data
par0 = load("parameters.mat").par0; %load initial parameters
sample_data = load("sample_data.mat").sample_data; %load the 'true' sample data 

%% set initial parameter values
num_params = [1,2,3,5,6,7,8,9,10,11,12,13,14,15,18,19,21,22,24,25]; %set trainable parameters
no_params = length(num_params); %number of trainable parameters
initial_par = par0(num_params); %get initial values of trainable parameters 

%% parametric bootstrap

noise = 0.1; %define amount of noise/error

%get 'true' data values
glucose = sample_data.glucose;
insulin = sample_data.insulin;
NEFA = sample_data.NEFA;
TG = sample_data.TG;

N = 50; %SET: number of new samples to make 
all_glucose = zeros(N, length(glucose));
all_insulin = zeros(N, length(insulin));
all_NEFA = zeros(N, length(NEFA));
all_TG = zeros(N, length(TG));
%resample data using a normal distribution
for i = 1:N
    all_glucose(i,:) = glucose + noise*randn(size(glucose)).*glucose;
    all_insulin(i,:) = insulin + noise*randn(size(insulin)).*insulin;
    all_NEFA(i,:) = NEFA + noise*randn(size(NEFA)).*NEFA;
    all_TG(i,:) = TG + noise*randn(size(TG)).*TG;
end

%make one big struct of all data samples
BW = sample_data.BW;
meal = sample_data.meal;

all_data(1,1) = sample_data; %the first struct must be the 'true' data 
for i = 2:N+1
    all_data(1,i).time_G = sample_data.time_G;
    all_data(1,i).time_I = sample_data.time_I;
    all_data(1,i).time_NEFA = sample_data.time_NEFA;
    all_data(1,i).time_TG = sample_data.time_TG;
    all_data(1,i).glucose = all_glucose(i-1,:);
    all_data(1,i).insulin = all_insulin(i-1,:);
    all_data(1,i).NEFA = all_NEFA(i-1,:);
    all_data(1,i).TG = all_TG(i-1,:);
    all_data(1,i).BW = BW;
    all_data(1,i).meal = meal;
end

save("bootstrap_data_test.mat", "all_data") %save sample data

%% fit model to different samples

time = sample_data.time_G;
lb = zeros(1,no_params); %set lower bounds to zero
ub = []; %set upper bounds 

all_par_opt = zeros(N,no_params);
for i = 1:N+1
    try
        result = Fit_M3al_Model(initial_par, all_data(i), time, lb, ub); %fit the model
        all_par_opt(i,:) = result.p_opt; %set the optimal parameters values in all_par_opt
    catch 
        continue
    end
end

save("bootstrap_par_test", "all_par_opt") %save all optimal parameters