%% This code is for conductivity of conductors affected by pellet microstructure
% By Howard Tu, 11/24/20
clc;clear;

%% Parameters
Good_LPS = 0.142;                    % The ionic conductivity of Good-LPS, in unit mS/cm
Bad_LPS  = 0.00675;                  % The ionic conductivity of Bad-LPS, in unit mS/cm
Good_LYC = 0.261;                    % The ionic conductivity of Good-LYC, in unit mS/cm
LLZO     = 0.5;                      % The ionic conductivity of LLZO, in unit mS/cm
sgm_le = [Good_LPS, Good_LYC]';      % The ionic conductivity of liquid electrolyte, in unit mS/cm
% sgm_le = Good_LPS;                    % The ionic conductivity of Good-LPS, in unit mS/cm
sgm_se = LLZO;                        % The ionic conductivity of the Oxide, in unit mS/cm
% The experimentally measured composite conductivity: 
% First column: volume fraction; Second column: conductivity
sgmExp_GoodLPS_LLZO = [0,Good_LPS; 0.15,0.212;   0.30,0.122; 0.45,0.0761;  0.6,0.0122]; 
% sgmExp_BadLPS_LLZO  = [0,0.0053;   0.15,Bad_LPS; 0.30,0.01;  0.45,0.00432; 0.6,0.00333]; 
sgmExp_GoodLYC_LLZO = [0,Good_LYC; 0.1,0.4; 0.15,0.372;  0.20,0.367; 0.3,0.26]; 

volFr  = linspace(0.01, 0.99, 99);              % The volume fraction of the particles
K_p    =  sgm_se ./ sgm_le;                 % The conductivity ratio between the SE and the LE

%% Solve the Bruggeman asymmetric equation for K_c:
sgm_cmp = zeros(length(sgm_le),length(volFr));          % The ionic conductivity of the composite, in unit mS/cm
syms y
for ifr = 1 : length(volFr)
    frac = volFr(ifr);
    K_c = zeros(length(sgm_le),1);              % The conductivity ratio between the composite and the LE
    for iLE = 1 : length(sgm_le)
        eqn          = y^3 - (1-frac)*(1-K_p(iLE))*y - K_p(iLE);
        ysol         = double(vpasolve( eqn == 0, y))';
        if real(ysol(1)) > 0
            K_c(iLE)   = ysol(1) .^ 3;
        else
            K_c(iLE)   = ysol(3) .^ 3;
        end
        
    end
    sgm_cmp(:,ifr) = K_c .* sgm_le;                    % The ionic conductivity of the composite, in unit mS/cm
end

%% Plot
figure(1)
for i = 1 : length(sgm_le)   
    plot(volFr,sgm_cmp(i,:));
    hold on
end
plot(sgmExp_GoodLPS_LLZO(:,1),sgmExp_GoodLPS_LLZO(:,2),'o')
hold on
plot(sgmExp_GoodLYC_LLZO(:,1),sgmExp_GoodLYC_LLZO(:,2),'o')
hold off
% legend('092520A','092420A','112720','092420B','102920A');
xlabel('volume fraction f');
ylabel('\sigma_c');

% figure(2)
% for i = 1 : length(sgm_le)   
%     plot(volFr,gradient(sgm_cmp(i,:),0.01));
%     hold on
% end
% hold off
% % legend('092520A','092420A','112720','092420B','102920A');
% xlabel('volume fraction');
% ylabel('dsgm_c / df');

% %% Define functions
% function myFun = mysol(x, frac, Kp)
%    myFun = x - (1-frac)*(1-Kp)*x^1/3 - Kp;
% end