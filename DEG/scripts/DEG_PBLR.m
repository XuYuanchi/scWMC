



filename = "/home/grads/fuzhowang2/suyanchi/scWMC/DEG/data/DEG.raw.csv";
iniData = readtable(filename, 'Delimiter', ',', 'ReadRowNames', false, 'ReadVariableNames', false);
minGenes = 0; minCells = 0; libararyflag = 0; logNormalize = 1;
proData = preprocessing(iniData, minCells, minGenes, libararyflag,logNormalize);
M = proData.data;
    
id = gene_selection(M);
    
M0 = M(id,:);
K = []; % the cluster numbers can be given by user.
numCores = 10; % defined by user
system_used = 'Mac';
accelerate = 0;
label = [];
[group,coph] = clusteing(iniData,M0,K,numCores,system_used,accelerate,label);
    
boundary_function = 3;
imputation_all = 1;
accelate = true;
PBLR_samp = PBLR_main(M,id,group,boundary_function,imputation_all,numCores,accelate);
filename = "/home/grads/fuzhowang2/suyanchi/scWMC/DEG/Results/PBLR.mat";
save(filename, 'PBLR_samp');
