%% scWMC DEMO: 
%% Clear all variables
clear;
addpath('lib/PROPACK','utils');
%% Load the data
% There are three datasets in the .mat file. There are the true data set,
% Drop-out data set, and the bulk data set.
load('data/demo_data.mat')
data_sc = data_dropout;
%% Parameters
Par.lam  = 0.8;
Par.rho  = 0.8  ;
Par.mu1  = 0.00001;
Par.mu2  = 0.00001;
Par.iter = 100; 
%% Run scWMC
dataRecovered = impute(data_sc, Par);
dataRecovered = max(dataRecovered, 0);
index         = find(data_sc);
dataRecovered(index) = data_sc(index); 
%% SAVE result
% writematrix(dataRecovered,'data/scWMC.csv');

Fro_error = norm(dataRecovered - data_true, 'fro');
disp(['****** L_2 error is ' num2str(Fro_error) '******']);
R    = corrcoef(dataRecovered, data_true);
disp(['******Pearson is ' num2str(R(1,2)) '******']);

%% Plot the results
gcf = figure(1);
set(gcf, 'Position', [100, 500, 1200, 300])
subplot(1,3,1)
imagesc(log10(data_true+1))
title('True Data')
axis off
subplot(1,3,2)
imagesc(log10(data_sc+1))
title('Drop-out Data')
axis off
subplot(1,3,3)
imagesc(log10(dataRecovered+1))
title('Imputed Data by scWMC')
axis off