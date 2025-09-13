function out = Q5_meal_model_parameters(p_opt,input_data)
%incorporate parameters to be estimated (p_opt) and fixed parameters in one 
%parameter vector for the Mixed Meal Model (M3al Model)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p_opt           - vector of parameter values to be optimised
% input_data      - struct of measured data needed to sepcify basal values
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for further information contact Shauna O'Donovan at
% s.d.odonovan@tue.nl
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%glucose + insulin parameters (EDES Rozendaal 2018) (fixed parameters are
%fixed to values used in Rozendaal et al 2018)

out(1)  = p_opt(1); %0.0105 %k1 rate constant for glucose stomach emptying (fast)[1/min]
out(2)  = p_opt(2); %0.28 %k2 rate constant for glucose appearence from gut [1/min]
out(3)  = p_opt(3); %6.07e-3 %k3 rate constant for suppresstion of hepatic glucose release by change of plasma glucose
out(4)  = 2.35e-4;  %k4 rate constant for suppression of hepatic glucose release by delayed (remote) insulin
out(5)  = p_opt(4); %0.0424;   %k5 rate constant for delayed insulin depedent uptake of glucose
out(6)  = p_opt(5); %2.2975;   %k6 rate constant for stimulation of insulin production by the change of plasma glucose concentration (beta cell funtion)
out(7)  = p_opt(6); %1.15;     %k7 rate constant for integral of glucose on insulin production (beta cell function)
out(8)  = p_opt(7); %7.27;     %k8 rate constant for the simulation of insulin production by the rate of change in plasma glucose concentration (beta cell function)
out(9)  = p_opt(8); %3.83e-2;  %k9 rate constant for outflow of insulin from plasma to interstitial space
out(10) = p_opt(9); %2.84e-1;  %k10 rate constant for degredation of insulin in remote compartment
out(11) = p_opt(10); %1.4;      %sigma shape factor (appearance of meal)
out(12) = p_opt(11); %13.2;     %Km michaelis-menten coefficient for glucose uptake
out(13) = p_opt(12); %input_data.glucose(1);%G_b basal plasma glucose [mmol/l]
out(14) = p_opt(13); %input_data.insulin(1); %I_PL/_b basal plasma glucose [microU/ml]
out(15) = p_opt(14); %0.043;    %basal hepatic glucose release

%triglyceride + NEFA parameters (new)
out(16) = 30;        %f_spill - fractional spillover of LPL derived NEFA
out(17) = 0.00045;   %k11 - rate coeficient LPL lipolysis sips/jelic models
out(18) = p_opt(15); %0.215;    %ATL_max maximum rate of ATL lipolysis in adipose tissues
out(19) = p_opt(16); %0.0385;   %K_ATL michealis menten coeficient for ATL lipolysis of store TG in adipose tissue
out(20) = 0.0713;    %k12 rate constanst for uptake of NEFA into tissues (currently insulin indenpendent)
out(21) = p_opt(17); %208;      %tau_LPL time delay for insulin stimulation of LPL lipolysis
out(22) = p_opt(18); %0.0088;   %k13 - rate constant for stomach emptying TG(very slow)
out(23) = 0.0163;    %k14 - rate constant for rate of TG appearance from gut
out(24) = p_opt(19); %1e-5;     %k15 coefficient for inhibition of TG secretion from liver by insulin
out(25) = p_opt(20); %0.0119; %k16 basal secretion of TG from liver

