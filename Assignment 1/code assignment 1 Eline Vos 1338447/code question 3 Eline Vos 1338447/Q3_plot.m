%this file makes the plots which can be found in section 3 of the report 

%%
%this part of the code was run to calculate the mean KS statistics for the
%dummies at lower amounts of ensembles with the implemented function. 
%Number of ensembles above 1000 were run separately with Q3_MPSA_no_sp.m 
%because these took quite long. 

% nr_ensembles_pos = [50, 100, 200, 300, 500, 800, 1000];
% means_dummies_low = zeros(1,length(nr_ensembles_pos));
% for i = 1:length(nr_ensembles)
%     nr_ensemble = nr_ensembles(i);
%     [Np, KS, mean_d] = Q3_MPSA_without_sp(5, nr_ensemble,1);
%     means_dummies(i) = mean_d;
% end
% save('mean_dummies_low.mat', "means_dummies_low") 

%%
%this part of the code was run to calculate the mean KS statistics for the
%dummies at ensembles where the scaling was substracted from the inital
%parameter value.

% nr_ensembles_neg = [100, 200, 300, 500, 800, 1000, 1500, 2000,2500];
% means_dummies_neg = zeros(1,length(nr_ensembles_neg));
% for i = 1:length(nr_ensembles)
%     nr_ensemble = nr_ensembles(i);
%     [Np, KS, mean_d] = Q3_MPSA_without_sp(5, nr_ensemble, 0);
%     means_dummies(i) = mean_d;
% end
% save('mean_dummies_neg.mat', "means_dummies_neg") 

%% plot with mean sensitivity dummies for varying amount of sample sets.
mean_dummy_low = load("mean_dummies_low.mat");
mean_dummy_1300 = load("mean_dummies_1300.mat");
mean_dummy_1600 = load("mean_dummies_1600.mat");
mean_dummy_2000 = load("mean_dummies_2000.mat");

dlow = mean_dummy_low.means_dummies;
d1300 = mean_dummy_1300.mean_dummies;
d1600 = mean_dummy_1600.mean_dummies;
d2000 = mean_dummy_2000.mean_dummies;

nr_ensembles = [50, 100, 200, 300, 500, 800, 1000, 1300, 1600, 2000];
d_all = [dlow, d1300, d1600, d2000];

figure(1)
plot(nr_ensembles, d_all, nr_ensembles, d_all, 'r*')
xlabel('number of ensembles')
ylabel('mean KS statistic')
%title('mean KS statistic for dummies for different amount of ensembles')

%%
mean_dummy_neg = load("mean_dummies_neg.mat");
dneg = mean_dummy_neg.means_dummies;

figure(2)
nr_ensembles_neg = [100, 200, 300, 500, 800, 1000, 1500, 2000, 2500];
plot(nr_ensembles_neg, dneg, nr_ensembles_neg, dneg, 'r*')
xlabel('number of ensembles')
ylabel('mean KS statistic')
%title('mean KS statistic for dummies for different amount of ensembles')


%% do the multiparametric sensitivity analyis for number of ensembles set to 1500
% [Np, KS, mean_d] = Q3_MPSA_without_sp(5, 1500, 1);
% [Np_sp, KS_sp, mean_d_sp] = Q3_MPSA_with_sp(5, 1500, 1);
% [Np_neg, KS_neg, mean_d_neg] = Q3_MPSA_without_sp(5, 1500, 0);
% [Np_nsp, KS_nsp, mean_d_nsp] = Q3_MPSA_with_sp(5, 1500, 0);

%save the results so you can look at the KS statiscs for each parameter
%more easily and you do not need to run above three lines again. 
%save('results_without_sp.mat', "Np", "KS", "mean_d") 
%save('results_with_sp.mat', "Np_sp", "KS_sp", "mean_d_sp")
%save('results_negative.mat', "Np_neg", "KS_neg", "mean_d_neg")
%save('results_negative_sp.mat', "Np_nsp", "KS_nsp", "mean_d_nsp")

%% plot KS statistics for model parameters and dummies with number of ensembles is 1500
results_without_sp = load("results_without_sp.mat");
results_with_sp = load("results_with_sp.mat");
results_neg = load("results_negative.mat");
results_nsp = load("results_negative_sp.mat");
nd = 5;

%plot without sample person characteristics 
Np = results_without_sp.Np;
KS = results_without_sp.KS;
X_param = 1:Np;
X_dummies = 1:nd;
mean_dummies = mean(KS(Np+1:end));

figure(3)
tiledlayout(1,2)

ax1 = nexttile;
scatter(X_param, KS(1:Np))
yline(mean_dummies)
ylabel('sensititivity (KS-statistic)', 'FontSize',12)
xlabel('parameter #', 'FontSize',12)

ax2 = nexttile;
scatter(X_dummies, KS(Np+1:end))
yline(mean_dummies)
xlabel('dummy #', 'FontSize',12)

linkaxes([ax1, ax2], 'y')

%plot with sample person characteristics)
Np_sp = results_with_sp.Np_sp;
KS_sp = results_with_sp.KS_sp;
X_param_sp = 1:Np_sp;
X_dummies_sp = 1:nd;
mean_dummies_sp = mean(KS_sp(Np_sp+1:end));

figure(4)
tiledlayout(1,2)

ax1 = nexttile;
scatter(X_param_sp, KS_sp(1:Np_sp))
yline(mean_dummies_sp)
ylabel('sensititivity (KS-statistic)', 'FontSize',12)
xlabel('parameter #', 'FontSize',12)

ax2 = nexttile;
scatter(X_dummies_sp, KS_sp(Np_sp+1:end))
yline(mean_dummies_sp)
xlabel('dummy #', 'FontSize',12)

linkaxes([ax1, ax2], 'y')

%plot without sample person characteristics and substraction
Np_neg = results_neg.Np_neg;
KS_neg = results_neg.KS_neg;
X_param_neg = 1:Np_neg;
X_dummies_neg = 1:nd;
mean_dummies_neg = mean(KS_neg(Np_neg+1:end));

figure(5)
tiledlayout(1,2)

ax1 = nexttile;
scatter(X_param_neg, KS_neg(1:Np_neg))
yline(mean_dummies_neg)
ylabel('sensititivity (KS-statistic)', 'FontSize',12)
xlabel('parameter #', 'FontSize',12)

ax2 = nexttile;
scatter(X_dummies_neg, KS_neg(Np_neg+1:end))
yline(mean_dummies_neg)
xlabel('dummy #', 'FontSize',12)

linkaxes([ax1, ax2], 'y')

%plot with sample person characteristics and substraction
Np_nsp = results_nsp.Np_nsp;
KS_nsp = results_nsp.KS_nsp;
X_param_nsp = 1:Np_nsp;
X_dummies_nsp = 1:nd;
mean_dummies_nsp = mean(KS_nsp(Np_nsp+1:end));

figure(6)
tiledlayout(1,2)

ax1 = nexttile;
scatter(X_param_nsp, KS_nsp(1:Np_nsp))
yline(mean_dummies_nsp)
ylabel('sensititivity (KS-statistic)', 'FontSize',12)
xlabel('parameter #', 'FontSize',12)

ax2 = nexttile;
scatter(X_dummies_nsp, KS_nsp(Np_nsp+1:end))
yline(mean_dummies_nsp)
xlabel('dummy #', 'FontSize',12)

linkaxes([ax1, ax2], 'y')


