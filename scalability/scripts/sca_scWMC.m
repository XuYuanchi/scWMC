%% scWMC DEMO: 
%% Clear all variables
clear;
addpath('code/lib/PROPACK','code/utils');
%% Load the data
% There are three datasets in the .mat file. There are the true data set,
% Drop-out data set, and the bulk data set.
%%
% 0.5593, 0.4153, 0.5830, 0.6067, 0.5935, 0.6875, 0.6312 
%%

    filename = "data/100k.csv";
    data_dropout = readtable(filename, 'Delimiter', ',', 'ReadRowNames', true, 'ReadVariableNames', true);
    data_dropout = table2array(data_dropout);
    data_dropout = process(data_dropout);
    data_dropout = single(data_dropout');
    [U, Z, V] = svd(data_dropout, "econ");
    
    %% dropout rate
    % disp(1-nnz(data_sc)/numel(data_sc));
    %% Parameters
    Par.lam  = 0.5;
    Par.rho  = 0.5  ;
    Par.mu1  = 0.00001;
    Par.mu2  = 0.00001;
    Par.iter = 5; 
    %% Run scWMC
    dataRecovered = impute(Z*V', Par);
    dataRecovered = (U*dataRecovered)';
    filename = "Results/scWMC.mat";
    save(filename, 'dataRecovered');

