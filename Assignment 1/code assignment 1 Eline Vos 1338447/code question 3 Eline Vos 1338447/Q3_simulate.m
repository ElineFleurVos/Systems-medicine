function feature = Simulate(phenotypic_data, parameters)
%input: struct with values of sample person characteristics and double with
%parameter values 

%output: feature, which is just a list with the output values of the four
%state variables in plasma at all timepoints

%%
%define time frame to look at
time = 0:1:500;
%define intial values and model constants needed for simulation of M3al Model model
[initial_values,constants] = M3al_Model_Initial(phenotypic_data,parameters);

%define global parameters for simulation
global t_saved G_PL_saved;
%initialise gloabl parameters
t_saved = 0;
G_PL_saved = phenotypic_data.glucose(1);

%specify options for ODE solver (Integrator function)
ODE_options = odeset('RelTol',1e-5,'OutputFcn',@integratorfunG);

%simulate model
[T,X] = ode45(@M3al_Model_ODE,time,initial_values,ODE_options,parameters,constants,phenotypic_data);

G_plasma = X(:,2);
I_plasma = X(:,4);
NEFA_plasma = X(:,9);
TG_plasma = X(:,13);

feature = [G_plasma; I_plasma; NEFA_plasma; TG_plasma];
