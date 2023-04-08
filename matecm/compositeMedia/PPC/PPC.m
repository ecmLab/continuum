%% This code is for conductivity measurement with Powder-Polymer composite electrolyte methode
% By Howard Tu, 06/14/21
clc;clear;

%% Parameters
% 1. Load the conductivity as a function of temperature
LLZO       = load('LLZO_conductivity.csv'); % 1st column:1000/T, unit 1/K; 2st column: log(sgm*T), unit S/cm * K
PEO       = load('PEO_conductivity.csv');  % 1st column:1000/T, unit 1/K; 2st column: sgm, unit S/cm
PEO_LLZO  = load('PEO_LLZO_conductivity.csv');  % 1st column:1000/T, unit 1/K; 2st column: sgm, unit S/cm
% 2. Convert all data unit to: 1st column:1000/T, unit 1/K; 2st column: log10(sgm), unit mS/cm
LLZO(:,2) = LLZO(:,2) + log10(LLZO(:,1));
PEO(:,2)  = log10(PEO(:,2)) + 3;
PEO_LLZO(:,2)  = log10(PEO_LLZO(:,2)) + 3;
log_sgm_se  = -1.6 * (PEO(:,1)-LLZO(1,1)) + LLZO(1,2);   % Get interplated LLZO value
% 3. Assign conductivities to SE and LE
sgm_se  = 10 .^ log_sgm_se;
sgm_le  = 10 .^ PEO(:,2);
% 4. volume fraction to be considered
volFr  = [0.2/1.2, 0.5/1.5];     % The volume fraction of the particles
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
plot(PEO(:,1),log_sgm_se,'-o',PEO(:,1),PEO(:,2),'-*',PEO_LLZO(:,1),PEO_LLZO(:,2),'-x')
legend('pure LLZO','pure PEO','wt% PEO:LLZO=1:0.5');
xlabel('1000/T [K^{-1}]');
ylabel('log(\sigma_c) [mS/cm]');

figure(2)
for i = 1 : length(volFr)   
    plot(PEO(:,1),log10(sgm_cmp(:,i)));
    hold on
end
plot(PEO(:,1),PEO(:,2),'-*',PEO_LLZO(:,1),PEO_LLZO(:,2),'-x')
hold off
xlabel('1000/T [K^{-1}]');
ylabel('log(\sigma_c) [mS/cm]');
