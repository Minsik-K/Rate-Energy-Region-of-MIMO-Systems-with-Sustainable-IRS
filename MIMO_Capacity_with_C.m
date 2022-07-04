function [Q, rate] = MIMO_Capacity_with_C(H, C, P)

M=size(H,2);

C = (C+C')/2;
lambda_min = max(real(eigs(C,1)), 0);  % lambda >= largest eigenvalue of C
lambda_max = (lambda_min+1)*10;
lambda_Max = lambda_max;

lambda_eps = 10^-7;

% while lambda_max - lambda_min > lambda_eps
for i = 1 : 1000
    lambda = (lambda_min + lambda_max) / 2;
    
    if lambda_Max-lambda < 2*lambda_eps
        disp('increase lambda_max');
    end
        
    TT=inv(sqrtm(lambda*eye(M)-C));
    [~, S, V]=svd(H*TT,'econ');
    Lambda = max(1-(1./diag(S)).^2,0);
    Q=TT*V*diag(Lambda)*V'*TT;
    
    if trace(Q)>P
        lambda_min = lambda;
    else
        lambda_max = lambda;
    end
    
    if lambda_max - lambda_min < lambda_eps
        break;
    end
end

rate = log2(det(eye(size(H,1)) + H*Q*H'));