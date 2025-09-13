%this file calculates for the different results which parameter
%sensitivities are above a certain set threshold, in this case 1 standard
%deviation above the mean of the dummie senstivities.

results_without_sp = load("results_without_sp.mat");
results_with_sp = load("results_with_sp.mat");
results_neg = load("results_negative.mat");
results_neg_sp = load("results_negative_sp.mat");

%analysis without sample person characteristics 
Np = results_without_sp.Np;
KS = results_without_sp.KS;
dummy_values = (KS(Np+1:end));
mean_dummies = mean(dummy_values); %calculate mean senitivity dummies
std_dummies = std(dummy_values); %calculate standard deviation sensitivity dummies

j = 1;
KS_params = (KS(1:Np));
for i = 1:length(KS_params)
    if KS_params(i) >= mean_dummies + 1*std_dummies %set threshold one standard deviation above mean
        params_sens(j) = i;
        j = j+1;
    end
end

%analysis with sample person characteristics 
Np_sp = results_with_sp.Np_sp;
KS_sp = results_with_sp.KS_sp;
dummy_values_sp = (KS_sp(Np_sp+1:end));
mean_dummies_sp = mean(dummy_values_sp); %calculate mean sensitivity dummies
std_dummies_sp = std(dummy_values_sp); %calculate standard deviation sensitivity dummies

j = 1;
KS_params_sp = (KS_sp(1:Np_sp));
for i = 1:length(KS_params_sp)
    if KS_params_sp(i) >= mean_dummies_sp + 1*std_dummies_sp %set threshold one standard deviation above mean
        params_sens_sp(j) = i;
        j = j+1;
    end
end

%analysis with parameters values below initial values without sample
%characteristics varied
Np_neg = results_neg.Np_neg;
KS_neg = results_neg.KS_neg;
dummy_values_neg = (KS_neg(Np_neg+1:end)); %calculate mean sensitivity dummies
std_dummies_neg = std(dummy_values_neg); %calculate standard deviation sensitivity dummies

mean_dummies_neg = mean(KS_neg(Np_neg+1:end));

j = 1;
KS_params_neg = (KS_neg(1:Np_neg));
for i = 1:length(KS_params_neg)
    if KS_params_neg(i) >= mean_dummies_neg + 1*std_dummies_neg %set threshold one standard deviation above mean
        params_sens_neg(j) = i;
        j = j+1;
    end
end

%analysis with parameters values below initial values with sample
%characteristics varied

Np_nsp= results_neg_sp.Np_nsp;
KS_nsp = results_neg_sp.KS_nsp;
dummy_values_nsp = (KS_nsp(Np_nsp+1:end)); %calculate mean sensitivity dummies
std_dummies_nsp = std(dummy_values_nsp); %calculate standard deviation sensitivity dummies

mean_dummies_nsp = mean(KS_nsp(Np_nsp+1:end));

j = 1;
KS_params_nsp = (KS_nsp(1:Np_nsp));
for i = 1:length(KS_params_nsp)
    if KS_params_nsp(i) >= mean_dummies_nsp + 1*std_dummies_nsp %set threshold one standard deviation above mean
        params_sens_nsp(j) = i;
        j = j+1;
    end
end

