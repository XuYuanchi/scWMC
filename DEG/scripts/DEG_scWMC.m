clear;
addpath('../../code/lib/PROPACK','../../code/utils');
%% load raw data
data_sc              = readtable('D:/MyWorkWorld/Imputation/scWMC/DEG/data/DEG.raw.csv','ReadRowNames',false,'ReadVariableNames',false);
% data_sc              = readtable('D:/MyWorkWorld/Imputation/scWMC/DEG/data/norm.csv','ReadRowNames',false,'ReadVariableNames',false);
data_sc              = table2array(data_sc);
%data_sc = process(data_sc);
% disp(1-nnz(data_sc)/numel(data_sc));
%% Parameters
Par.lam  = 0.4;
Par.rho  = 0.4;
Par.mu1  = 0.000001;
Par.mu2  = 0.000001;
Par.iter = 10; 
%% imputing
dataRecovered = impute(data_sc, Par);
dataRecovered = max(dataRecovered, 0);
index         = find(data_sc);
dataRecovered(index) = data_sc(index);
%% save
filename = "D:/MyWorkWorld/Imputation/scWMC/DEG/Results/scWMC.mat";
save(filename, 'dataRecovered');
