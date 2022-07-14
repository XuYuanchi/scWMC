% addpath(genpath('./'));
% addpath PROPACK;
% addpath Rank;

data_name = ["Deng", "Petropoulos"];
for i=2:2
    filename = "D:/MyWorkWorld/Imputation/scWMC/ICT/data/" + data_name(1,i) + "_filter.csv";
    iniData = readtable(filename, 'Delimiter', ',', 'ReadRowNames', true, 'ReadVariableNames', true);
    minGenes = 0; minCells = 0; libararyflag = 0; logNormalize = 1;
    proData = preprocessing(iniData, minCells, minGenes, libararyflag,logNormalize);
    M = proData.data;
    
    id = gene_selection(M);
    
    M0 = M(id,:);
    K = []; % the cluster numbers can be given by user.
    numCores = 1; % defined by user
    system_used = 'Mac';
    accelerate = 0;
    label = [];
    [group,coph] = clusteing(iniData,M0,K,numCores,system_used,accelerate,label);
    
    boundary_function = 3;
    imputation_all = 1;
    accelate = true;
    PBLR_samp = PBLR_main(M,id,group,boundary_function,imputation_all,numCores,accelate);
    filename = "D:/MyWorkWorld/Imputation/scWMC/ICT/Results/" + data_name(1,i) + "/PBLR.mat";
    save(filename, 'PBLR_samp');
end