function [X] = alm_WMC(M, Omega, Q_U, Q_V, maxIter,mu1,mu2,isDisp)

%mu1 = 0.01;mu2 = 0.01;
tol = 1e-3;tolProj = 1e-5;
[m, n] = size(M);

Omega_c = ones(m,n) - Omega;

% I = eye(n);


% initialize

X = zeros( m, n);
E = zeros( m, n);
Z = zeros( m, n);
Y1 = zeros( m, n);
Y2 = zeros( m, n);
R1 = zeros( m, n);
R2 = zeros( m, n);
% invM = inv(mu1*I+mu2*Q*Q');

invQ_U = inv(Q_U);
invD = inv(Q_U'*Q_U);
mu   = mu1/mu2;
A    = mu*invD;
B    = Q_V*Q_V';

% dnorm = norm(M, 'fro');
% tolProj = tolProj * dnorm;
total_svd = 0;


iter = 0;
converged = false;
% stopCriterion = 1;
sv = 5;
svp = sv;
max_sv = min(m,n);
while ~converged       
    iter = iter + 1;
    % solve the primal problem by alternative projection
    primal_converged = false;
    primal_iter = 0;
    sv = min(sv + round(max_sv * 0.1),max_sv);
    while primal_converged == false && primal_iter <10
        
        C      = mu*invD*(M-E) + invD*Y1/mu2 + invQ_U*(Z+Y2)*Q_V';
        
        temp_X = sylvester(A, B, C);
        
%         temp_X = (mu1*M+mu1*E-Y1+mu2*Z*Q'+Y2*Q')*invM;
        E = Omega_c.*(M-temp_X+Y1/mu1);
        
        
        
        if choosvd(n, sv) == 1
            [U, Sig, V] = lansvd(Q_U*temp_X*Q_V - (1/mu2)*Y2, sv, 'L');
        else
            [U, Sig, V] = svd(Q_U*temp_X*Q_V - (1/mu2)*Y2, 'econ');
        end
        diagSig = diag(Sig);
        svp = length(find(diagSig > 1/mu2));
        if svp < sv
            sv = min(svp + 1, max_sv);
        else
            sv = min(svp + round(0.05*max_sv), max_sv);
        end
        Z = U(:,1:svp)*diag(diagSig(1:svp)-1/mu2)*V(:,1:svp)';    
        
        
%         if norm(X - temp_X, 'fro') < tolProj*norm(X,'fro')
        error_inner = norm(log10(X + 1) - log10(temp_X + 1), 'fro')/(m*n);
        if error_inner < 1e-4
            primal_converged = true;
        end
        X = temp_X;
        primal_iter = primal_iter + 1;
        total_svd = total_svd + 1;
               
    end
     
    temp_R1 = M - X - E;
    temp_R2 = Z-Q_U*X*Q_V;        
    Y1 = Y1 + mu1*temp_R1;
    Y2 = Y2 + mu2*temp_R2;
    
    %% stop Criterion    
%     stopCriterion1 = norm(temp_R1 - R1, 'fro');
%     stopCriterion2 = norm(temp_R2 - R2, 'fro');
    stopCriterion1 = norm(log10(temp_R1+1) - log10(R1+1), 'fro');
    stopCriterion2 = norm(log10(temp_R2+1) - log10(R2+1), 'fro');
    if stopCriterion1 < tol && stopCriterion2 < tol
        converged = true;
    end    
    R1 = temp_R1;
    R2 = temp_R2;
    
   if isDisp
       disp(['Iteration' num2str(iter) ' #svd ' num2str(total_svd) ' Rank(X) ' num2str(svp)...
            ' stopCriterion1 ' num2str(stopCriterion1)  ' stopCriterion2 ' num2str(stopCriterion2)]);
   end
    
    if ~converged && iter >= maxIter
        %disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end


