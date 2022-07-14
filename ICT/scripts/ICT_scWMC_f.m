%% scWMC DEMO: 
%% Clear all variables
clear;
addpath('../../code/lib/PROPACK','../../code/utils');

data_name = ["Deng", "Petropoulos"];
% 0.396
% 0.3565


i=2;
%    filename = "D:/MyWorkWorld/Imputation/scWMC/ICT/data/" + data_name(i) + "_filter.csv"; 
   filename = "D:/MyWorkWorld/Imputation/scWMC/ICT/data/" + data_name(i) + "_filter.csv";
   
   data_dropout = readtable(filename, 'Delimiter', ',', 'ReadRowNames', true, 'ReadVariableNames', true);
   
   data_sc = table2array(data_dropout);
   % data_sc(all(data_sc  == 0,2), :) = [];
   data_sc = process(data_sc);
   
   %% fater
   [U, Z, V] = svd(data_sc, "econ");
   
   %% dropout rate
   disp(1-nnz(data_sc)/numel(data_sc));
   %% Parameters
    Par.lam  = 0.4;
    Par.rho  = 0.4;
    Par.mu1  = 0.001;
    Par.mu2  = 0.001;
    Par.iter = 100;
   %% Run scWMC
    dataRecovered = impute(Z*V', Par);
    dataRecovered = (U*dataRecovered)';
%     dataRecovered = max(dataRecovered, 0);
    index         = find(data_sc);
    dataRecovered(index) = data_sc(index);
    filename = "D:/MyWorkWorld/Imputation/scWMC/ICT/Results/"+ data_name(i) + "/scWMC.mat";
%     filename = "D:/MyWorkWorld/Imputation/scWMC/ICT/Results/Petropoulos/scWMC.mat";
    save(filename, 'dataRecovered');


