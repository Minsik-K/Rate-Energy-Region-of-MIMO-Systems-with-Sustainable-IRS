function [rate, energy, Q, alpha, obj_traj] = Proposed_Algorithm(H, R, T, P, eta, w, Mod, convergence)

Iter = 100;
eps=10^-4;

Nt = size(H,2);  % no. of transmit antennas
Nr = size(H,1);  % no. of receive antennas
M = size(R,2);  % no. of IRS shifts

alpha = ones(M,1);
Q = MIMO_Capacity_with_C(H, zeros(Nt,Nt), P);
D = w*eta*T'*Q*T;

alpha =  CD( H, R, T, Q, alpha, D, Mod);
C = w*eta*T*diag(1-abs(alpha))*T';
H_tilde = H + R*diag(alpha)*T';
   
old_obj=0;
for iter=1:Iter
    
    % Optimizing Q based on [11]
    Q = MIMO_Capacity_with_C(H_tilde, C, P);
    D = w*eta*T'*Q*T;
    
    % Optimizing PHI while fixing Q
    [alpha, ~, ~, ~] = CD( H, R, T, Q, alpha, D, Mod);
    
    H_tilde = H + R*diag(alpha)*T';
    C = w*eta*T*diag(1-abs(alpha))*T';
    
    % Objective function = Rate + w*Harvested Energy
    rate = real(log2(det(eye(Nr)+H_tilde*Q*H_tilde')));
    weighted_energy = real(trace(C*Q));

    obj_traj(iter) = rate + weighted_energy;
    
    % Stop criterion
    if abs(obj_traj(iter)-old_obj)<eps && convergence == 0
        break
    else
        old_obj=obj_traj(iter);
    end
    
end

rate = real(log2(det(eye(Nr)+H_tilde*Q*H_tilde')));
if w==0 % if w is zero, all IRS elements are in active mode
    energy = 0;
else
    energy = real(trace(C*Q))/w;
end



function [alpha, obj, rate, weighted_energy] = CD(H, R, T, Q, alpha, D, Mod)

Iter = 100;

Nt = size(H,2);  % no. of transmit antennas
Nr = size(H,1);  % no. of receive antennas
M = size(R,2);  % no. of IRS shifts

eps = 10^-7;

H_tilde = H + R*diag(alpha)*T';

old_obj = 0;
for iter=1:Iter

    for m=1:M
        H_tilde = H_tilde - alpha(m)*R(:,m)*T(:,m)';

        Am = eye(Nr) + H_tilde*Q*H_tilde' + R(:,m)*T(:,m)'*Q*T(:,m)*R(:,m)';
        Bm = R(:,m)*T(:,m)'*Q*H_tilde';
        
        temp=trace(inv(Am)*Bm)'+10^-10; % To avoid division by zero
        temp=temp/abs(temp);

        alpha(m) = temp;
        
        % Projection to PSK constellation
        if Mod ~= 0
            alpha(m) = pskmod(pskdemod(alpha(m),Mod,pi/Mod),Mod,pi/Mod);
        end

        varphi0 = real(log2(det(eye(Nr) + H_tilde*Q*H_tilde')));
        varphi1 = real(log2(det(Am + alpha(m)*Bm + alpha(m)*Bm')));
        
        if varphi1 < real(D(m,m)) + varphi0
            alpha(m)=0;
        end
    
        H_tilde = H_tilde + alpha(m)*R(:,m)*T(:,m)';
    end    
        
    rate = real(log2(det(eye(Nr)+H_tilde*Q*H_tilde')));
    weighted_energy = real(trace(D*diag(1-abs(alpha))));
    
    obj = rate + weighted_energy;
    
    if obj-old_obj<eps 
        break
    else
        old_obj=obj;
    end            
    
end
