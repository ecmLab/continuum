%% This code is for conductivity of conductors affected by pellet microstructure
% By Howard Tu, 11/24/20
clc;clear;

%% Parameters
sgm_le =  0.05*10 .^ linspace(0,2,100)';      % The ionic conductivity of liquid electrolyte, in unit S/cm
sgm_se =  0.5;                          % The ionic conductivity of solid electrolyte, in unit S/cm
volFr  =  [0.2,0.47];              % The volume fraction of the particles

K_p    =  sgm_se ./ sgm_le;                 % The conductivity ratio between the SE and the LE

%% Solve the Bruggeman asymmetric equation for K_c:
%(K_c - K_p)/(K_c^1/3*(1-K_p)) = 1 - f
% Method 1
% K_c1  = zeros(length(sgm_le),1);              % The conductivity ratio between the composite and the LE
% for iLE = 1 : length(sgm_le)
%     K_c1(iLE)     = fsolve(@(x) mysol(x, frac, K_p(iLE)), 0.001);
% end
% Method 2
sgm_cmp = zeros(length(sgm_le),length(volFr));          % The ionic conductivity of the composite, in unit S/cm
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
    sgm_cmp(:,ifr) = K_c .* sgm_le;                    % The ionic conductivity of the composite, in unit S/cm
end

%% Plot
figure(1)
plot(sgm_le,sgm_cmp(:,1),'-b');
hold on
plot(sgm_le,sgm_cmp(:,2),'-r');
hold on
% plot(log10(sgm_le),log10(sgm_cmp(:,3)),'-g');
% hold on
plot(sgm_le,sgm_le,'--k');
hold off
% figure(2)
% plot(K_p,K_c,'--x');
%
% %% Define functions
% function myFun = mysol(x, frac, Kp)
%    myFun = x - (1-frac)*(1-Kp)*x^1/3 - Kp;
% end