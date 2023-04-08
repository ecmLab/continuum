%% This code is for conductivity of conductors affected by pellet microstructure
% Used for comparing the paper "A high conductivity oxide?sulfide composite lithium superionic conductor"
% By Howard Tu, 11/24/20
clc;clear;

%% Parameters
% sgm_le = [0.003765, 0.005908, 0.01066, 0.01339, 0.07996]';    
sgm_lps  = 0.16;                           % The ionic conductivity of LPS, in unit mS/cm
sgm_llzo = 0.4;                            % The ionic conductivity of LLZO, in unit mS/cm
den_lps  = 1.83;                           % The density of LPS, in unit g/cm^3
den_llzo = 5.11;                           % The density of LLZO, in unit g/cm^3

% The experimentally measured value
% 1st column: weight fraction of LLZO; 2nd column: relative density; 3th: composite conductivity;
expV = [0.0,96.2,0.1412; 0.1,90.1,0.2938; 0.2,83.6,0.3640; 0.25,80.3,0.4052; 0.3,78.4,0.5336; 0.35,77.2,0.3955; 0.4,75.8,0.3277; 0.6,74.4,0.2163];
% volFr  = linspace(0.01, 0.99, 99);            % The volume fraction of the LLZO particles
K_p     =  sgm_llzo ./ sgm_lps;                 % The conductivity ratio between the SE and the LE

%% Solve the Bruggeman asymmetric equation for K_c when zero porosity:
volFr         =  zeros(length(expV),1);         % The volume fraction of the LLZO particles at zero porosity
volFr(2:end)  =  1/(1+(1-expV(2:end,1))/expV(2:end,1)*den_llzo/den_lps); 
sgm_cmp = zeros(length(sgm_lps),length(volFr));          % The ionic conductivity of the composite, in unit mS/cm
syms y
for ifr = 1 : length(volFr)
    frac = volFr(ifr);
    K_c = zeros(length(sgm_lps),1);              % The conductivity ratio between the composite and the LE
    for iLE = 1 : length(sgm_lps)
        eqn          = y^3 - (1-frac)*(1-K_p(iLE))*y - K_p(iLE);
        ysol         = double(vpasolve( eqn == 0, y))';
        if real(ysol(1)) > 0
            K_c(iLE)   = ysol(1) .^ 3;
        else
            K_c(iLE)   = ysol(3) .^ 3;
        end
        
    end
    sgm_cmp(:,ifr) = K_c .* sgm_lps;                    % The ionic conductivity of the composite, in unit mS/cm
end

%% Solve the Bruggeman equation with given porosity:


%% Plot
figure(1)
for i = 1 : length(sgm_lps)   
    plot(volFr,sgm_cmp(i,:));
    hold on
end
plot(expValue(:,1),expValue(:,2),'o')
hold off
% legend('092520A','092420A','112720','092420B','102920A');
xlabel('volume fraction');
ylabel('sgm_c');

figure(2)
for i = 1 : length(sgm_lps)   
    plot(volFr,gradient(sgm_cmp(i,:),0.01));
    hold on
end
hold off
% legend('092520A','092420A','112720','092420B','102920A');
xlabel('volume fraction');
ylabel('dsgm_c / df');

% %% Define functions
% function myFun = mysol(x, frac, Kp)
%    myFun = x - (1-frac)*(1-Kp)*x^1/3 - Kp;
% end