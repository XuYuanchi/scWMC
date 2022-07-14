function [M_rec] = impute(Ori_P, Par)

%% calculate column and row space and its orthogonal complement
[U , ~, V]  = svd(Ori_P,'econ');
[m, n]      = size(Ori_P);

%% calculate orthogonal projections
P_U      = U * U';
P_U_     = eye(m) - P_U;
P_V      = V * V';
P_V_     = eye(n) - P_V;

%% calculate projection matrix
lam      = Par.lam;
rho      = Par.rho;
Q_U      = lam * P_U + P_U_;
Q_V      = rho * P_V + P_V_;

%% complete the gene matrix by WMC
[m, n]   = size(Ori_P);
Omega    = zeros(m,n);
Omega(Omega~= Ori_P) = 1; 
iter_num = Par.iter;
mu1      = Par.mu1;
mu2      = Par.mu2;

M_rec    = alm_WMC(Ori_P, Omega, Q_U, Q_V, iter_num, mu1, mu2, 1);
end