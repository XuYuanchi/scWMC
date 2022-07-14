% addpath(genpath('./'));
% addpath PROPACK;
% addpath Rank;

data_name = ["sc_CELseq2", "sc_10x", "sc_Droseq", "sc_10x_5cl", "sc_Celseq2_5cl_p1", "sc_Celseq2_5cl_p2", "sc_Celseq2_5cl_p3"];
for j =6:10
for i=3:6
    filename = "D:/MyWorkWorld/Imputation/scWMC/clustering/data/" + data_name(i) + ".csv";
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
    filename = "D:/MyWorkWorld/Imputation/scWMC/clustering/PBLR/" + num2str(j) + "/" + data_name(i);
    save(filename, 'PBLR_samp');
end
end