%This file creates the three mock datasets that will be analysed in
%question 7

noise = 0.1; %amount of noise in mock data (10%)
%Define the new sampling schedules
ts0 = [1,60,130,240,360,480];
ts1 = [1,50,100,150,200,240];
ts2 = [1,100,200,400];
ts3 = [1,100,200,300];

%% Specify phenotypic traits and meal 
%specify phenotypic traits
sample_person.glucose = 5;    %fasting glucose (mmol/l)
sample_person.insulin = 18;   %fasting insulin (uIU/ml)
sample_person.TG      = 1.3;  %fasting plasma triglyceride (mmol/l)
sample_person.NEFA    = 0.33; %fasting plasma NEFA(mmol/l)
sample_person.BW      = 84.2; %body weight (kg)

%specify composition of meal being simulated
sample_person.meal.G  = 75000; %mass of glucose in meal (mg)
sample_person.meal.TG = 60000; %mass of triglyceride in meal (mg)

%load inital parameters
parameters = load("parameters.mat").par0;

%% simulate the initial model
time = 0:1:500; %time span for model simulation
[T, X] = Q2_simulate_meal_model(sample_person,parameters,time);

%get the inital simulated data for the four state variables
G_data = X(:,2)';
I_data = X(:,4)';
NEFA_data = X(:,9)';
TG_data = X(:,13)';

sample_data = load("sample_data.mat").sample_data;

%general data, which is the same for each sampling schedulde
BW = sample_person.BW;
meal = sample_person.meal;

%% create sample data set for given sample schedule, but then wih noise

mock_data(1,1).time_G = ts0;
mock_data(1,1).time_I = ts0;
mock_data(1,1).time_TG = ts0;
mock_data(1,1).time_NEFA = ts0;
%Explanation what is done below: noise is added on top of the true value. The random normalized value is
%thus both multiplied with the true value and the noise level as the addition
%is done with respect to the true data.
mock_data(1,1).glucose = G_data(ts0) + G_data(ts0).*noise.*randn(size(G_data(ts0)));
mock_data(1,1).insulin = I_data(ts0) + I_data(ts0).*noise.*randn(size(I_data(ts0)));
mock_data(1,1).NEFA = NEFA_data(ts0) + NEFA_data(ts0).*noise.*randn(size(NEFA_data(ts0)));
mock_data(1,1).TG = TG_data(ts0) + TG_data(ts0).*noise.*randn(size(TG_data(ts0)));

%% create sample data set for sample schedule 1

mock_data(1,2).time_G = ts1;
mock_data(1,2).time_I = ts1;
mock_data(1,2).time_TG = ts1;
mock_data(1,2).time_NEFA = ts1;
%Explanation what is done below: noise is added on top of the true value. The random normalized value is
%thus both multiplied with the true value and the noise level as the addition
%is done with respect to the true data.
mock_data(1,2).glucose = G_data(ts1) + G_data(ts1).*noise.*randn(size(G_data(ts1)));
mock_data(1,2).insulin = I_data(ts1) + I_data(ts1).*noise.*randn(size(I_data(ts1)));
mock_data(1,2).NEFA = NEFA_data(ts1) + NEFA_data(ts1).*noise.*randn(size(NEFA_data(ts1)));
mock_data(1,2).TG = TG_data(ts1) + TG_data(ts1).*noise.*randn(size(TG_data(ts1)));

%% create sample data set for sample schedule 2

mock_data(1,3).time_G = ts2;
mock_data(1,3).time_I = ts2;
mock_data(1,3).time_TG = ts2;
mock_data(1,3).time_NEFA = ts2;
mock_data(1,3).glucose = G_data(ts2) + G_data(ts2).*noise.*randn(size(G_data(ts2)));
mock_data(1,3).insulin = I_data(ts2) + I_data(ts2).*noise.*randn(size(I_data(ts2)));
mock_data(1,3).NEFA = NEFA_data(ts2) + NEFA_data(ts2).*noise.*randn(size(NEFA_data(ts2)));
mock_data(1,3).TG = TG_data(ts2) + TG_data(ts2).*noise.*randn(size(TG_data(ts2)));

%% create sample data set for sample schedule 3

mock_data(1,4).time_G = ts3;
mock_data(1,4).time_I = ts3;
mock_data(1,4).time_TG = ts3;
mock_data(1,4).time_NEFA = ts3;
mock_data(1,4).glucose = G_data(ts3) + G_data(ts3).*noise.*randn(size(G_data(ts3)));
mock_data(1,4).insulin = I_data(ts3) + I_data(ts3).*noise.*randn(size(I_data(ts3)));
mock_data(1,4).NEFA = NEFA_data(ts3) + NEFA_data(ts3).*noise.*randn(size(NEFA_data(ts3)));
mock_data(1,4).TG = TG_data(ts3) + TG_data(ts3).*noise.*randn(size(TG_data(ts3)));

%% add bodyweight and meal to all at ones
for i = 1:4
    mock_data(1,i).BW = BW;
    mock_data(1,i).meal = meal;
end

save("mock_data_test", "mock_data")

