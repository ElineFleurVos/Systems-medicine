function errors = Q7_simulate_meal_model(phenotypic_data,parameters,time)
% Simulate M3al Model for a given parameter set. Also plot the sample data 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phenotypic_data - struct of information specifying fasting glucose,
%                   insulin, TG, and NEFA. Must also sepcify body weight (BW)
%                   and meal composition. 
% parameters      - parameter values to be simulated.
% time            - timespan for ode simulation

%% simulate model for given parameter set

%define intial values and model constants needed for simulation of M3al Model model
[initial_values,constants] = M3al_Model_Initial(phenotypic_data,parameters);

%define global parameters for simulation
global t_saved G_PL_saved;
%initialise gloabl parameters
t_saved = 0;
G_PL_saved = phenotypic_data.glucose(1);

%specify options for ODE solver (Integrator function)
ODE_options = odeset('RelTol',1e-5,'OutputFcn',@integratorfunG);

par_true = load("parameters.mat").par0;
par0 = parameters(1,:);
par1 = parameters(2,:);
par2 = parameters(3,:);
par3 = parameters(4,:);

%simulate model
[T_true, X_true] = ode45(@M3al_Model_ODE,time,initial_values,ODE_options,par_true,constants,phenotypic_data);
[T0,X0] = ode45(@M3al_Model_ODE,time,initial_values,ODE_options,par0,constants,phenotypic_data);
[T1,X1] = ode45(@M3al_Model_ODE,time,initial_values,ODE_options,par1,constants,phenotypic_data);
[T2,X2] = ode45(@M3al_Model_ODE,time,initial_values,ODE_options,par2,constants,phenotypic_data);
[T3,X3] = ode45(@M3al_Model_ODE,time,initial_values,ODE_options,par3,constants,phenotypic_data);

%calculates SSE for different sampling schedules
X_true_all = [X_true(:,2);X_true(:,4);X_true(:,9);X_true(:,13)];
X0_all = [X0(:,2);X0(:,4);X0(:,9);X0(:,13)];
X1_all = [X1(:,2);X1(:,4);X1(:,9);X1(:,13)];
X2_all = [X2(:,2);X2(:,4);X2(:,9);X2(:,13)];
X3_all = [X3(:,2);X3(:,4);X3(:,9);X3(:,13)];
error1 = sum((X_true_all'-X0_all').^2);
error2 = sum((X_true_all'-X1_all').^2);
error3 = sum((X_true_all'-X2_all').^2);
error4 = sum((X_true_all'-X3_all').^2);
errors = [error1, error2, error3, error4];

%% Generate figure
 plot_col1 = [0, 0.4470, 0.7410]; 
 plot_col2 = [0.8500 0.3250 0.0980];
 plot_col3 = [0.9290 0.6940 0.1250];
 plot_col4 = [0.4940 0.1840 0.5560];
 plot_col5 = [0,0,0];

 mock_data = load("mock_data.mat").mock_data;
 ss0 = mock_data(1);
 ss1 = mock_data(2);
 ss2 = mock_data(3);
 ss3 = mock_data(4);
 
 subplot(2,2,1)
 %plot plasma glucose
 %plot the model fits for the tree sample schedules
 plot(T0,X0(:,2),'Color', plot_col1,'LineWidth',1.5);
 hold on
 plot(T1,X1(:,2),'Color', plot_col2,'LineWidth',1.5);
 hold on
 plot(T2,X2(:,2),'Color', plot_col3,'LineWidth',1.5);
 hold on
 plot(T3,X3(:,2),'Color', plot_col4,'LineWidth',1.5);
 hold on
 plot(T_true,X_true(:,2),'Color', plot_col5,'LineWidth',1.5);
 hold on
 %also plot the sample data 
 scatter(ss0.time_G, ss0.glucose,30, plot_col1, 'filled')
 hold on
 scatter(ss1.time_G, ss1.glucose,30, plot_col2, 'filled')
 hold on 
 scatter(ss2.time_G, ss2.glucose,30, plot_col3, 'filled')
 hold on 
 scatter(ss3.time_G, ss3.glucose,30, plot_col4, 'filled')
 ylabel('glucose (mmol/l)');
 title('plasma Glucose')
 xlabel('time (mins)')
 legend('ss0','ss1', 'ss2', 'ss3', 'true model')
 
 subplot(2,2,2)
 %plasma insulin
 %plot the model fits for the tree sample schedules
 plot(T0,X0(:,4),'Color',plot_col1,'LineWidth',1.5);
 hold on
 plot(T1,X1(:,4),'Color',plot_col2,'LineWidth',1.5);
 hold on
 plot(T2,X2(:,4),'Color', plot_col3,'LineWidth',1.5);
 hold on
 plot(T3,X3(:,4),'Color', plot_col4,'LineWidth',1.5);
 hold on
 plot(T_true,X_true(:,4),'Color', plot_col5,'LineWidth',1.5);
 hold on
 %also plot the sample data 
 scatter(ss0.time_I, ss0.insulin,30, plot_col1, 'filled')
 hold on 
 scatter(ss1.time_I, ss1.insulin,30, plot_col2, 'filled')
 hold on 
 scatter(ss2.time_I, ss2.insulin,30, plot_col3, 'filled')
 hold on 
 scatter(ss3.time_I, ss3.insulin,30, plot_col4, 'filled')
 ylabel('insulin (uIU/ml)');
 title('plasma insulin')
 legend('ss0','ss1', 'ss2', 'ss3', 'true model')
 
 subplot(2,2,3)
 %plasma triglyceride
 %plot the model fits for the tree sample schedules
 plot(T0,X0(:,13),'Color',plot_col1,'LineWidth',1.5);
 hold on
 plot(T1,X1(:,13),'Color',plot_col2,'LineWidth',1.5);
 hold on
 plot(T2,X2(:,13),'Color', plot_col3,'LineWidth',1.5);
 hold on
 plot(T3,X3(:,13),'Color', plot_col4,'LineWidth',1.5);
 hold on
 plot(T_true,X_true(:,13),'Color', plot_col5,'LineWidth',1.5);
 hold on
 %also plot the sample data 
 scatter(ss0.time_TG, ss0.TG,30, plot_col1, 'filled')
 hold on 
 scatter(ss1.time_TG, ss1.TG,30, plot_col2, 'filled')
 hold on 
 scatter(ss2.time_TG, ss2.TG,30, plot_col3, 'filled')
 hold on 
 scatter(ss3.time_TG, ss3.TG,30, plot_col4, 'filled')
 ylabel('TG (mmol/l)');
 title('plasma TG')
 xlabel('time (mins)')
 legend('ss0','ss1', 'ss2', 'ss3', 'true model')
 
 subplot(2,2,4)
 %plasma NEFA
 %plot the model fits for the tree sample schedules
 plot(T0,X0(:,9),'Color',plot_col1,'LineWidth',1.5);
 hold on
 plot(T1,X1(:,9),'Color',plot_col2,'LineWidth',1.5)
 hold on
 plot(T2,X2(:,9),'Color', plot_col3,'LineWidth',1.5);
 hold on
 plot(T3,X3(:,9),'Color', plot_col4,'LineWidth',1.5);
 hold on
 plot(T_true,X_true(:,9),'Color', plot_col5,'LineWidth',1.5);
 hold on
 %also plot the sample data 
 scatter(ss0.time_NEFA, ss0.NEFA,30, plot_col1, 'filled')
 hold on
 scatter(ss1.time_NEFA, ss1.NEFA,30, plot_col2, 'filled')
 hold on 
 scatter(ss2.time_NEFA, ss2.NEFA,30, plot_col3, 'filled')
 hold on 
 scatter(ss3.time_NEFA, ss3.NEFA,30, plot_col4, 'filled')
 ylabel('NEFA (mmol/l)');
 title('plasma NEFA')
 xlabel('time (mins)')
 legend('ss0','ss1', 'ss2', 'ss3', 'true model')
 