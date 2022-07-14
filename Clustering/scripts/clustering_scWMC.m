%% scWMC DEMO: 
%% Clear all variables
clear;
addpath('../code/lib/PROPACK','../code/utils');
%% Load the data
% There are three datasets in the .mat file. There are the true data set,
% Drop-out data set, and the bulk data set.
data_name = ["sc_CELseq2", "sc_10x", "sc_Droseq", "sc_10x_5cl", "sc_Celseq2_5cl_p1", "sc_Celseq2_5cl_p2", "sc_Celseq2_5cl_p3"];
%%
% 0.5593, 0.4153, 0.5830, 0.6067, 0.5935, 0.6875, 0.6312 
%%
for j =1:10
for i=3:6
    filename = "D:/MyWorkWorld/Imputation/scWMC/Clustering/data/" + data_name(i) + ".mat";
%     data_dropout = readtable(filename, 'Delimiter', ',', 'ReadRowNames', true, 'ReadVariableNames', true);
    load(filename);
    label = double(categorical(label));
    cluster_num = length(unique(label));
    data_sc = process(data_sc);
%     data_sc = table2array(data_dropout);
    
    %% dropout rate
    % disp(1-nnz(data_sc)/numel(data_sc));
    %% Parameters
    Par.lam  = 0.6;
    Par.rho  = 0.6;
    Par.mu1  = 0.001;
    Par.mu2  = 0.001;
    Par.iter = 50; 
    %% Run scWMC
    dataRecovered = impute(data_sc, Par);
    dataRecovered = max(dataRecovered, 0);
    index         = find(data_sc);
    dataRecovered(index) = data_sc(index);
    idx = kmeans(dataRecovered.', cluster_num);
    NMI = nmi(label, idx);
    disp(NMI);
    filename = "D:/MyWorkWorld/Imputation/scWMC/Clustering/scWMC/" + num2str(j) + "/" + data_name(i) + "_mm.mat";
    save(filename, 'dataRecovered');
end
end