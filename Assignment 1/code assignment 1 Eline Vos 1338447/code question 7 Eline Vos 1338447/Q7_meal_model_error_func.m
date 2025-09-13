function out = Q7_meal_model_error_func(p_opt,input_data,time)
%Cost function to fit Mixed Meal Model to measured data. The error between 
%the Mixed Meal Model (M3al Model) simulation for a given parameter set and 
%supplied measured meal challenge test data is calculated. Additional
%constraints are applied to ensure certain physioloigical constraints are
%met; 
%i)   All glucose contained in the meal appears in 4 hours following the meal. 
%ii)  All triglyceride contained in the meal appears in 12 hours following the meal.
%iii) At fasting (t=0 mins)NEFA uptake into tissues is in steady state (equal to measured value)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p_opt        - vector of parameter values for which error is being
%                calculated.
% input_data   - struct of measured challenge test data for caculating error
%              - mean and standard deviation values are required for
%                glucose and insulin (need to extend to TG and NEFA)
%              - vector of sampling time points are also required.
%              - any additional variables required for model simulation.
% time         - time span for model simulation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for further information contact Shauna O'Donovan at
% s.d.odonovan@tue.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simulate model for given parameter set

%form full parameter vector for simulation
parameters = Q7_meal_model_parameters(p_opt,input_data);
%define intial values and model constants needed for simulation of eDES model
[initial_values,constants]=M3al_Model_Initial(input_data,parameters);

%define global parameters for simulation
global t_saved G_PL_saved;
%initialise gloabl parameters
t_saved = 0;
G_PL_saved = input_data.glucose(1);

%specify options for ODE solver (Integrator function)
ODE_options = odeset('RelTol',1e-5,'OutputFcn',@integratorfunG);

%simulate model
[T,X] = ode45(@M3al_Model_ODE,time,initial_values,ODE_options,parameters,constants,input_data);

%% Calculate error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model fit error - data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
measured_time_G=ismember(T,input_data.time_G);
G_err = (X(measured_time_G,2)' - input_data.glucose);

%%insulin error
measured_time_I=ismember(T,input_data.time_I);
I_err = (X(measured_time_I,4)' - input_data.insulin);

%NEFA error
measured_time_NEFA=ismember(T,input_data.time_NEFA);
NEFA_err = 10*(X(measured_time_NEFA,9)' - input_data.NEFA); %multiply the NEFA errors with 10

%TG error
measured_time_TG=ismember(T,input_data.time_TG);
TG_err = 10*(X(measured_time_TG,13)' - input_data.TG); %multiply the TG errors with 10

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% total error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out=[G_err,I_err,NEFA_err,TG_err];