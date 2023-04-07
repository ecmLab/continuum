% parameter
clc;clear;

%% Parameters
% Constants
FF    = 96485.33;              % Faraday constant, in unit s*A/mol
RT    = 8.314*300;             % RT constant, in unit J/mol
Vm    = 12.97;                 % Molar volume of Lithium, in unit cm^3/mol
F_RT  = FF/RT;                 % Constant: F/RT, unit 1/V
Vm_F  = Vm/FF * 3600 *1000;    % Constant: Molar volume of Li / Faraday constant, in unit um^3/(Hour*nA)
ifg   = 0;                     % Figure plot flag
% Electrochemistry parameters
i_exc = 1.3;  % The exchange current density, mA/cm^2
alfa  = 0.5;  % Reaction rate
% Geometry of the model
lx    =  1e-3;                      % length of the model, in unit um
lyz   =  1e-3;                      % height of the model, in unit um

%% Load data
% Meaning of each column in the csv file:
% 1st: time, in unit second;                                            2st: defGrad_zz: total deformation gradient along z direction;
% 3st: effPlastic_strain: effective plastic strain;                     4st: eff_strain: effective total strain;
% 5st: fp_zz: plastic deformation gradient along z direction;           6st: hard: the harding strength, in unit Pa;
% 7st: logStrain_zz: logarithmic elastic strain along z direction;      8st: mandel_zz: Mandel stress along z direction, in unit Pa;
% 9st: matl_ts_min: ;                                                   10st: strain_rate: plastic strain rate;
% 11st: strain_zz: total strain along z direction;                      12st: stress_zz: Cauchy stress along z direction, in unit Pa;
% 13st: u_z: average displacement at the surface of z=zmax, in unit um; 14st: von_mises: Von-mises stress, in unit Pa;
% 15st: yield_strength: yield strength, in unit Pa;
% Tension test under constant strain rate
tension_strnR2000  = csvread('tension_constTrueStrainRateR2000.csv',1,0);
tension_strnR300   = csvread('tension_constTrueStrainRateR300.csv',1,0);
tension_strnR30    = csvread('tension_constTrueStrainRateR30.csv',1,0);
tension_strnR4     = csvread('tension_constTrueStrainRateR4.csv',1,0);
% compression test under constant strain rate
compression_strnR2000  = csvread('compression_constTrueStrainRateR2000.csv',1,0);
compression_strnR300   = csvread('compression_constTrueStrainRateR300.csv',1,0);
compression_strnR30    = csvread('compression_constTrueStrainRateR30.csv',1,0);
% tension test under constant stress
tension_strsC10  = csvread('tension_constTrueStressC10.csv',1,0);
tension_strsC8   = csvread('tension_constTrueStressC8.csv',1,0);
tension_strsC6   = csvread('tension_constTrueStressC6.csv',1,0);
% compression test under constant stress
compression_strsC10  = csvread('compression_constTrueStressC10.csv',1,0);
compression_strsC8   = csvread('compression_constTrueStressC8.csv',1,0);
compression_strsC6   = csvread('compression_constTrueStressC6.csv',1,0);

%% Data plot
% Plot true-stress vs. time along tension direction 
% stress in unit MPa, time in unit seconds
ifg = ifg + 1; 
figure(ifg);
plot(tension_strnR2000(:,1),tension_strnR2000(:,11)*100,'-o'); 
hold on
plot(tension_strnR300(:,1),tension_strnR300(:,11)*100,'-o'); 
% hold on
% plot(tension_strnR03(:,1),tension_strnR03(:,12)/10^6,'-o'); 
% hold on
% plot(tension_strnR004(:,1),tension_strnR004(:,12)/10^6,'-o'); 
hold off

% Plot true-stress vs. true-strain along tension and compression direction under constant strain rate condition
% Strain in unit %, stress in unit MPa
ifg = ifg + 1; 
figure(ifg);
% Tension
plot(tension_strnR2000(:,11)*100,tension_strnR2000(:,12)/10^6,'-ob'); 
hold on
plot(tension_strnR300(:,11)*100,tension_strnR300(:,12)/10^6,'-xr'); 
hold on
plot(tension_strnR30(:,11)*100,tension_strnR30(:,12)/10^6,'-ok'); 
hold on
plot(tension_strnR4(:,11)*100,tension_strnR4(:,12)/10^6,'-.'); 
% Compression
hold on
plot(compression_strnR2000(:,11)*100,compression_strnR2000(:,12)/10^6,'-ob'); 
hold on
plot(compression_strnR300(:,11)*100,compression_strnR300(:,12)/10^6,'-xr'); 
hold on
plot(compression_strnR30(:,11)*100,compression_strnR30(:,12)/10^6,'-ok'); 
hold off
axis([-26,26,-1.5,1.5]);
title('Constant Strain rate');

% Plot true-stress vs. true-strain along tension and compression direction under constant stress condition
% Strain in unit %, stress in unit MPa
ifg = ifg + 1; 
figure(ifg);
% Tension
plot(tension_strsC10(:,11)*100,tension_strsC10(:,12)/10^6,'-ob'); 
hold on
plot(tension_strsC8(:,11)*100,tension_strsC8(:,12)/10^6,'-xr'); 
hold on
plot(tension_strsC6(:,11)*100,tension_strsC6(:,12)/10^6,'-ok'); 
hold on
% Compression
plot(compression_strsC10(:,11)*100,compression_strsC10(:,12)/10^6,'-ob'); 
hold on
plot(compression_strsC8(:,11)*100,compression_strsC8(:,12)/10^6,'-xr'); 
hold on
plot(compression_strsC6(:,11)*100,compression_strsC6(:,12)/10^6,'-ok'); 
hold off
axis([-30,30,-1.5,1.5]);
title('Constant Stress');
