function [rate, energy, Q, alpha, obj] = full_search(H, R, T, P, eta, w ,Mod)

Nt = size(H,2);  % no. of transmit antennas
Nr = size(H,1);  % no. of receive antennas
M = size(R,2);  % no. of IRS shifts

set_digit = pskmod(0:(Mod)-1,Mod,pi/Mod);
set_digit = [0, set_digit];

rate_vec = zeros(1, 3^M);
weighted_energy_vec = zeros(1, 3^M);
obj_vec = zeros(1, 3^M);

for iter=1:(Mod+1)^M

    temp = dec2base(iter-1,Mod+1,M);
    for i = 1:M
        alpha(i,1) = set_digit(uint8(str2num(temp(i)))+1);
    end
    
    H_tilde = H + R*diag(alpha)*T';
    
    C = w*eta*T*diag(1-abs(alpha))*T';
    Q = MIMO_Capacity_with_C(H_tilde, C, P);
    
    rate_vec(iter) = real(log2(det(eye(Nr)+H_tilde*Q*H_tilde')));
    weighted_energy_vec(iter) = real(trace(C*Q));
    
    obj_vec(iter) = rate_vec(iter) + weighted_energy_vec(iter);
end

[obj, index] = max(obj_vec);
rate = rate_vec(index);
energy = weighted_energy_vec(index)/w;

temp = dec2base(index-1,Mod+1,M);
for i = 1:M
    alpha(i,1) = set_digit(uint8(str2num(temp(i)))+1);
end





