function [h_TR, h_IR, h_TI] = IRS_channel(Nt, Nr, M, option, d_x, beta_TR, beta_TI, beta_IR)
%
% Our channel model is based on following paper.
% M. Cui, G. Zhang, and R. Zhang, ``Secure wireless communication via intelligent reflecting 
% surface," \emph{IEEE Wireless Commun. Lett.}, vol. 8, no. 5, pp. 1410-1414, Oct. 2019.
%

if nargin<5
    option=0;  % Rician fiading LOS components
end    

zeta0 = 10^-3;

d0 = 1;
d_TR = 50;
d_vertical = 5;
d_TI = sqrt(d_x^2+d_vertical^2);
d_IR = sqrt((d_TR-d_x)^2+d_vertical^2);

sigma = sqrt(10^-8);  % noise variance -80 dBm

% the spatial correlation matrix
% R_c=zeros(Nt,Nt);
% for i=1:Nt
%     for j=1:Nt
%         R_c(i,j)=0.95^(abs(i-j));
%     end
% end 

% -------------------- Tx - Rx channel --------------------
K=0;
LOS = sqrt(K/(K+1));
NLOS = sqrt(1/(K+1));

if option==0
    c_L = ones(Nr,Nt);               % Rician fiading LOS components with same phase
else    
	c_L = exp(1i*rand(Nr,Nt)*2*pi);  % Rician fiading LOS components with random phase
end    
c_NL = (randn(Nr,Nt)+1i*randn(Nr,Nt))/sqrt(2);
h_TR = sqrt(zeta0*(d0/d_TR)^beta_TR) * (LOS*c_L + NLOS*c_NL) / sigma;


% -------------------- IRS channels ------------------------
K=3;
LOS = sqrt(K/(K+1));
NLOS = sqrt(1/(K+1));
% -------------------- Tx - IRS channel --------------------
%c_L = randn(1,Nt)+1i*randn(1,Nt); c_L=c_L./abs(c_L);
if option==0
    c_L = ones(Nt,M);
else    
	c_L = exp(1i*rand(Nt,M)*2*pi);
end    
c_NL = (randn(Nt,M)+1i*randn(Nt,M))/sqrt(2);
h_TI = sqrt(zeta0*(d0/d_TI)^beta_TI) * (LOS*c_L + NLOS*c_NL);


% -------------------- IRS - Rx channel --------------------
%c_L = randn(1,M)+1i*randn(1,M); c_L=c_L./abs(c_L);
if option==0
    c_L = ones(Nr,M);
else    
	c_L = exp(1i*rand(Nr,M)*2*pi);
end    
c_NL = (randn(Nr,M)+1i*randn(Nr,M))/sqrt(2);
h_IR = sqrt(zeta0*(d0/d_IR)^beta_IR) * (LOS*c_L + NLOS*c_NL) / sigma;

end

