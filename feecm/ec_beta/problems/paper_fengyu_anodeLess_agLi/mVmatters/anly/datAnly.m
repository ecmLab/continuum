clc; clear; myFigureSetting;

ifg = 0;
%Geometric parameter
xSE = 20;      % width of the SE, in unit um
ySE = 40;      % thickness of the SE, in unit um
dw  = 2;        % width of the defect located at the top middle, in unit um
dh  = 2;        % length of the defect, in unit um
%Electrochemical parameters
i_exc = 13;       % The exchange current density, based on the ASR=2 Ohm, in unit mA/cm^2
F_RT  = 0.03868;  % The combined constant F/RT when T=300K, in unit 1/mV
a0    = 0.5;      % Reaction rate for the charge transfer reaction, set as symmetric anodic and cathodic reaction for now
sgm   = [0.01, 0.1, 1, 10];  % The ionic conductivity, set as variable for parameter studies, in unit mS/cm

%% Read data from files
% Relation of filenames and variable names: "sgm001" and "sgm-2" means sgm = 10^-2 = 0.01 mS/cm 
% Load anode potential
andPot_sgm001 = csvread('../rst/andPot_sgm-2.csv',1,0);
andPot_sgm01  = csvread('../rst/andPot_sgm-1.csv',1,0);
andPot_sgm0   = csvread('../rst/andPot_sgm0.csv',1,0);
andPot_sgm1   = csvread('../rst/andPot_sgm1.csv',1,0);
% Load anode current
andCrnt_sgm001 = csvread('../rst/andCrnt_sgm-2.csv',1,0);
andCrnt_sgm01  = csvread('../rst/andCrnt_sgm-1.csv',1,0);
andCrnt_sgm0   = csvread('../rst/andCrnt_sgm0.csv',1,0);
andCrnt_sgm1   = csvread('../rst/andCrnt_sgm1.csv',1,0);

%% Data analysis
% Compute the tangent of bottom cosine curve: y_prime = -pi*dh/dw*sin(pi/dw*x)
np = length(andCrnt_sgm001(:,4));
xx = andCrnt_sgm001(:,4);  yp001 = zeros(np,1);    yp001(abs(xx)<dw/2) = -pi*dh/dw *sin(pi/dw*xx(abs(xx)<dw/2));
xx = andCrnt_sgm01(:,4);   yp01  = zeros(np,1);    yp01(abs(xx)<dw/2)  = -pi*dh/dw *sin(pi/dw*xx(abs(xx)<dw/2));
xx = andCrnt_sgm0(:,4);    yp0   = zeros(np,1);    yp0(abs(xx)<dw/2)   = -pi*dh/dw *sin(pi/dw*xx(abs(xx)<dw/2));
xx = andCrnt_sgm1(:,4);    yp1   = zeros(np,1);    yp1(abs(xx)<dw/2)   = -pi*dh/dw *sin(pi/dw*xx(abs(xx)<dw/2));
% Compute normal current density from numerical results 
 % i_n = ix*nx + iy*ny = ix*ty - iy*tx = (ix*y_prime - iy)/sqrt(1+y_prime^2)
 % Note: This is only for benckmark test to check the numerical stability and accuracy
 % Tagent vector: [nx, ny]; Normal vector: [tx, ty]
nrmCrnt_sgm001 = (andCrnt_sgm001(:,1).*yp001 - andCrnt_sgm001(:,2)) ./ sqrt(1+yp001.^2);
nrmCrnt_sgm01  = (andCrnt_sgm01(:,1).*yp01   - andCrnt_sgm01(:,2))  ./ sqrt(1+yp01.^2);
nrmCrnt_sgm0   = (andCrnt_sgm0(:,1).*yp0     - andCrnt_sgm0(:,2))   ./ sqrt(1+yp0.^2);
nrmCrnt_sgm1   = (andCrnt_sgm1(:,1).*yp1     - andCrnt_sgm1(:,2))   ./ sqrt(1+yp1.^2);
% Compute normal current density from the analytic Butler-Volmer relation:
 % i_n = i_exc * (exp(alpha*eta*F/RT)-exp(-alpha*eta*F/RT))
nrmCrnt0_sgm001 = i_exc * (exp(a0*F_RT*andPot_sgm001(:,2)) - exp(-a0*F_RT*andPot_sgm001(:,2)));
nrmCrnt0_sgm01  = i_exc * (exp(a0*F_RT*andPot_sgm01(:,2)) - exp(-a0*F_RT*andPot_sgm01(:,2)));
nrmCrnt0_sgm0   = i_exc * (exp(a0*F_RT*andPot_sgm0(:,2)) - exp(-a0*F_RT*andPot_sgm0(:,2)));
nrmCrnt0_sgm1   = i_exc * (exp(a0*F_RT*andPot_sgm1(:,2)) - exp(-a0*F_RT*andPot_sgm1(:,2)));

%% Plotting
% Plot the potential
ifg = ifg + 1;
figure(ifg)
plot(andPot_sgm001(:,3),andPot_sgm001(:,2),'-b',andPot_sgm01(:,3),andPot_sgm01(:,2),'-k',andPot_sgm0(:,3),andPot_sgm0(:,2),'-r');
legend('0.01mS/cm','0.1mS/cm','1mS/cm');
% Plot the current
ifg = ifg + 1;
figure(ifg)
plot(andCrnt_sgm001(:,4),nrmCrnt_sgm001,'-b',andCrnt_sgm01(:,4),nrmCrnt_sgm01,'-k',andCrnt_sgm0(:,4),nrmCrnt_sgm0,'-r');
hold on
plot(andCrnt_sgm001(:,4),nrmCrnt0_sgm001,'-b',andCrnt_sgm01(:,4),nrmCrnt0_sgm01,'-k',andCrnt_sgm0(:,4),nrmCrnt0_sgm0,'-r');
legend('0.01mS/cm','0.1mS/cm','1mS/cm');
hold off
