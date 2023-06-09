clc; clear; myFigureSetting;

ifg = 0;
%Geometric parameter
xSE = 40;      % width of the SE, in unit um
ySE = 40;      % thickness of the SE, in unit um
dw  = 2;        % width of the defect located at the top middle, in unit um
dh  = 2;        % length of the defect, in unit um
%Electrochemical parameters
F_RT  = 0.03868;  % The combined constant F/RT when T=300K, in unit 1/mV
a0    = 0.5;      % Reaction rate for the charge transfer reaction, set as symmetric anodic and cathodic reaction for now
% For LPS vs Li-metal system
i_exc = 13;       % The exchange current density, 13 is based on the ASR=2 Ohm, in unit mA/cm^2; LPS/Li metal interface from paper
sgmLPS=1;         % The ionic conductivity of LPS, set to be 1mS/cm
% For change variables
iSgm  = [0.01, 0.1, 1];  % The ionic conductivity, set as variable for parameter studies, in unit mS/cm
iExc  = [0.13, 1.3, 13]; % The exchange current density, set as variable for parameter studies, in unit mA/cm^2
iVlt  = [0,-1,-2,-3,-4,-5]; % The equilibrium chemical potential of the interface, in unit mV; 0 corresponds to Li metal anode; -5mV is the high Li content AgLi SS 
iAg   = [0,-1,-2,-3,-4,-5]; % The equilibrium voltage of the AgLi particle in the surface defect, in unit mV; 0 corresponds to Li metal anode;
iPrs  = [1 ,2, 3];  % The pressure in unit MPa, 1MPa corresponds to 0.175mV equilibrium potential change

%% Read data from files
% Relation of filenames and variable names: "sgm001" and "sgm-2" means sgm = 10^-2 = 0.01 mS/cm 
% Load anode potential
andPot_sgm001 = csvread('../rst/andPot_sgm-2.csv',1,0);
andPot_sgm01  = csvread('../rst/andPot_sgm-1.csv',1,0);
andPot_sgm0   = csvread('../rst/andPot_sgm0.csv',1,0);
andPot_sgm1   = csvread('../rst/andPot_sgm1.csv',1,0);
andPot_exc01  = csvread('../rst/andPot_exc-1.csv',1,0);
andPot_exc0   = csvread('../rst/andPot_exc0.csv',1,0);
andPot_exc1   = csvread('../rst/andPot_exc1.csv',1,0);
andPot_vlt0   = csvread('../rst/andPot_vlt0.csv',1,0);
andPot_vlt1  = csvread('../rst/andPot_vlt-1.csv',1,0);
andPot_vlt2   = csvread('../rst/andPot_vlt-2.csv',1,0);
andPot_vlt3   = csvread('../rst/andPot_vlt-3.csv',1,0);
andPot_vlt4   = csvread('../rst/andPot_vlt-4.csv',1,0);
andPot_vlt5   = csvread('../rst/andPot_vlt-5.csv',1,0);
andPot_Ag0    = csvread('../rst/andPot_Ag0.csv',1,0);
andPot_Ag1    = csvread('../rst/andPot_Ag-1.csv',1,0);
andPot_Ag2    = csvread('../rst/andPot_Ag-2.csv',1,0);
andPot_Ag3    = csvread('../rst/andPot_Ag-3.csv',1,0);
andPot_Ag4    = csvread('../rst/andPot_Ag-4.csv',1,0);
andPot_Ag5    = csvread('../rst/andPot_Ag-5.csv',1,0);

% Load anode current
andCrnt_sgm001 = csvread('../rst/andCrnt_sgm-2.csv',1,0);
andCrnt_sgm01  = csvread('../rst/andCrnt_sgm-1.csv',1,0);
andCrnt_sgm0   = csvread('../rst/andCrnt_sgm0.csv',1,0);
andCrnt_sgm1   = csvread('../rst/andCrnt_sgm1.csv',1,0);
andCrnt_exc01  = csvread('../rst/andCrnt_exc-1.csv',1,0);
andCrnt_exc0   = csvread('../rst/andCrnt_exc0.csv',1,0);
andCrnt_exc1   = csvread('../rst/andCrnt_exc1.csv',1,0);
andCrnt_vlt0   = csvread('../rst/andCrnt_vlt0.csv',1,0);
andCrnt_vlt1   = csvread('../rst/andCrnt_vlt-1.csv',1,0);
andCrnt_vlt2   = csvread('../rst/andCrnt_vlt-2.csv',1,0);
andCrnt_vlt3   = csvread('../rst/andCrnt_vlt-3.csv',1,0);
andCrnt_vlt4   = csvread('../rst/andCrnt_vlt-4.csv',1,0);
andCrnt_vlt5   = csvread('../rst/andCrnt_vlt-5.csv',1,0);
andCrnt_Ag0    = csvread('../rst/andCrnt_Ag0.csv',1,0);
andCrnt_Ag1    = csvread('../rst/andCrnt_Ag-1.csv',1,0);
andCrnt_Ag2    = csvread('../rst/andCrnt_Ag-2.csv',1,0);
andCrnt_Ag3    = csvread('../rst/andCrnt_Ag-3.csv',1,0);
andCrnt_Ag4    = csvread('../rst/andCrnt_Ag-4.csv',1,0);
andCrnt_Ag5    = csvread('../rst/andCrnt_Ag-5.csv',1,0);

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
nrmCrnt0_exc01  = iExc(1) * (exp(a0*F_RT*andPot_exc01(:,2)) - exp(-a0*F_RT*andPot_exc01(:,2)));
nrmCrnt0_exc0   = iExc(2) * (exp(a0*F_RT*andPot_exc0(:,2)) - exp(-a0*F_RT*andPot_exc0(:,2)));
nrmCrnt0_exc1   = iExc(3) * (exp(a0*F_RT*andPot_exc1(:,2)) - exp(-a0*F_RT*andPot_exc1(:,2)));
nrmCrnt0_vlt0   = i_exc * (exp(a0*F_RT*(andPot_vlt0(:,2) - iVlt(1))) - exp(-a0*F_RT*(andPot_vlt0(:,2) - iVlt(1))));
nrmCrnt0_vlt1   = i_exc * (exp(a0*F_RT*(andPot_vlt1(:,2) - iVlt(2))) - exp(-a0*F_RT*(andPot_vlt1(:,2) - iVlt(2))));
nrmCrnt0_vlt2   = i_exc * (exp(a0*F_RT*(andPot_vlt2(:,2) - iVlt(3))) - exp(-a0*F_RT*(andPot_vlt2(:,2) - iVlt(3))));
nrmCrnt0_vlt3   = i_exc * (exp(a0*F_RT*(andPot_vlt3(:,2) - iVlt(4))) - exp(-a0*F_RT*(andPot_vlt3(:,2) - iVlt(4))));
nrmCrnt0_vlt4   = i_exc * (exp(a0*F_RT*(andPot_vlt4(:,2) - iVlt(5))) - exp(-a0*F_RT*(andPot_vlt4(:,2) - iVlt(5))));
nrmCrnt0_vlt5   = i_exc * (exp(a0*F_RT*(andPot_vlt5(:,2) - iVlt(6))) - exp(-a0*F_RT*(andPot_vlt5(:,2) - iVlt(6))));
% For the AgLi particle case, things are a little bit complicated: zero voltage at the flat area, but AgLi voltage at the interface
eta0 = andPot_Ag0(:,2);      eta0(abs(andPot_Ag0(:,3))-dw/2<1) = andPot_Ag0(abs(andPot_Ag0(:,3))-dw/2<1,2) - iAg(1);
nrmCrnt0_Ag0    = i_exc * (exp(a0*F_RT*eta0) - exp(-a0*F_RT*eta0));
eta1 = andPot_Ag1(:,2);      eta1(abs(andPot_Ag1(:,3))-dw/2<1) = andPot_Ag1(abs(andPot_Ag1(:,3))-dw/2<1,2) - iAg(2);
nrmCrnt0_Ag1    = i_exc * (exp(a0*F_RT*eta1) - exp(-a0*F_RT*eta1));
eta2 = andPot_Ag2(:,2);      eta2(abs(andPot_Ag2(:,3))-dw/2<1) = andPot_Ag2(abs(andPot_Ag2(:,3))-dw/2<1,2) - iAg(3);
nrmCrnt0_Ag2    = i_exc * (exp(a0*F_RT*eta2) - exp(-a0*F_RT*eta2));
eta3 = andPot_Ag3(:,2);      eta3(abs(andPot_Ag3(:,3))-dw/2<1) = andPot_Ag3(abs(andPot_Ag3(:,3))-dw/2<1,2) - iAg(4);
nrmCrnt0_Ag3    = i_exc * (exp(a0*F_RT*eta3) - exp(-a0*F_RT*eta3));
eta4 = andPot_Ag4(:,2);      eta4(abs(andPot_Ag4(:,3))-dw/2<1) = andPot_Ag4(abs(andPot_Ag4(:,3))-dw/2<1,2) - iAg(5);
nrmCrnt0_Ag4    = i_exc * (exp(a0*F_RT*eta4) - exp(-a0*F_RT*eta4));
eta5 = andPot_Ag5(:,2);      eta5(abs(andPot_Ag5(:,3))-dw/2<1) = andPot_Ag5(abs(andPot_Ag5(:,3))-dw/2<1,2) - iAg(6);
nrmCrnt0_Ag5    = i_exc * (exp(a0*F_RT*eta5) - exp(-a0*F_RT*eta5));

%% Plotting
% Plot the potential vs normal current at sgm=0.01mS/cm
ifg = ifg + 1;
figure(ifg)
plot(andPot_sgm01(:,3),andPot_sgm01(:,2),'-b')
title('Voltage at SE/Anode Interface, in mV');
ifg = ifg + 1;
figure(ifg)
plot(andCrnt_sgm01(:,4),nrmCrnt0_sgm01,'-r')
title('Current at SE/Anode Interface, in mA/cm^2');

% Plot potential at ASR= Ohm*cm^2 and varying sgm=[0.01,0.1,1]mS/cm
ifg = ifg + 1;
figure(ifg)
plot(andPot_sgm01(:,3),andPot_sgm01(:,2),'-b',andPot_sgm0(:,3),andPot_sgm0(:,2),'-k',andPot_sgm1(:,3),andPot_sgm1(:,2),'-r');
legend('0.1mS/cm','1mS/cm-LPS','10mS/cm');
title('Voltage at SE/Anode Interface, in mV');
% Plot the current at ASR= 2 Ohm*cm^2 and varying sgm=[0.01,0.1,1]mS/cm
ifg = ifg + 1;
figure(ifg)
% plot(andCrnt_sgm001(:,4),nrmCrnt_sgm001,'-b',andCrnt_sgm01(:,4),nrmCrnt_sgm01,'-k',andCrnt_sgm0(:,4),nrmCrnt_sgm0,'-r');
% hold on
plot(andCrnt_sgm01(:,4),nrmCrnt0_sgm01,'-b',andCrnt_sgm0(:,4),nrmCrnt0_sgm0,'-k',andCrnt_sgm1(:,4),nrmCrnt0_sgm1,'-r');
legend('0.1mS/cm','1mS/cm-LPS','10mS/cm');
title('Current at SE/Anode Interface, in mA/cm^2');
% hold off

% Plot potential at sgm= 1 mS/cm and varying ASR= [0.13, 1.3, 13]Ohm*cm^2
ifg = ifg + 1;
figure(ifg)
plot(andPot_exc01(:,3),andPot_exc01(:,2),'-b',andPot_exc0(:,3),andPot_exc0(:,2),'-k',andPot_exc1(:,3),andPot_exc1(:,2),'-r');
legend('200 Ohm*cm^2','20 Ohm*cm^2','2 Ohm*cm^2-LPS');
title('Voltage at SE/Anode Interface, in mV');
% axis([-20,20,-1,-0.4]);
% Plot the current at sgm= 1 mS/cm and varying ASR= [0.13, 1.3, 13]Ohm*cm^2
ifg = ifg + 1;
figure(ifg)
plot(andCrnt_exc01(:,4),nrmCrnt0_exc01,'-b',andCrnt_exc0(:,4),nrmCrnt0_exc0,'-k',andCrnt_exc1(:,4),nrmCrnt0_exc1,'-r');
legend('200 Ohm*cm^2','20 Ohm*cm^2','2 Ohm*cm^2-LPS');
title('Current at SE/Anode Interface, in mA/cm^2');

% Plot potential at sgm= 1 mS/cm, ASR= 2 Ohm*cm^2 and varying anode equilibrium potential= [0, -1, -2, -3, -4 ,-5] mV
ifg = ifg + 1;
figure(ifg)
plot(andPot_vlt0(:,3),andPot_vlt0(:,2),'-b',  andPot_vlt1(:,3),andPot_vlt1(:,2),'-r', andPot_vlt2(:,3),andPot_vlt2(:,2),'-k');
legend('0 mV','1 mV','2 mV');
title('Voltage at SE/Anode Interface, in mV');
% Plot the current at sgm= 1 mS/cm, ASR= 2 Ohm*cm^2 and varying anode equilibrium potential= [0, -1, -2, -3, -4 ,-5] mV
ifg = ifg + 1;
figure(ifg)
plot(andCrnt_vlt0(:,4),nrmCrnt0_vlt0,'-b',  andCrnt_vlt1(:,4),nrmCrnt0_vlt1,'-r', andCrnt_vlt2(:,4),nrmCrnt0_vlt2,'-k');
legend('0 mV','1 mV','2 mV');
title('Current at SE/Anode Interface, in mA/cm^2');

% Plot potential at sgm= 1 mS/cm, ASR= 2 Ohm*cm^2 and varying AgLi particle equilibrium potential in the surface defect = [0, -1, -2, -3, -4 ,-5] mV
ifg = ifg + 1;
figure(ifg)
plot(andPot_Ag0(:,3),andPot_Ag0(:,2),'-b',andPot_Ag1(:,3),andPot_Ag1(:,2),'-k',andPot_Ag2(:,3),andPot_Ag2(:,2),'-r');
legend('0 mV','1 mV','2 mV');
title('Voltage at SE/Anode Interface, in mV');
% Plot the current at sgm= 1 mS/cm, ASR= 2 Ohm*cm^2 and varying AgLi particle equilibrium potential in the surface defect = [0, -1, -2, -3, -4 ,-5] mV
ifg = ifg + 1;
figure(ifg)
plot(andCrnt_Ag0(:,4),nrmCrnt0_Ag0,'-b',andCrnt_Ag1(:,4),nrmCrnt0_Ag1,'-k',andCrnt_Ag2(:,4),nrmCrnt0_Ag2,'-r');
legend('0 mV','1 mV','2 mV');
title('Current at SE/Anode Interface, in mA/cm^2');
