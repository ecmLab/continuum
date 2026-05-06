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
% Tension test under constant strain rate at temperature T = 198K
tension_strnT198R2000  = csvread('tension_constTrueStrainRateT198R2000.csv',1,0);
tension_strnT198R300   = csvread('tension_constTrueStrainRateT198R300.csv',1,0);
tension_strnT198R30    = csvread('tension_constTrueStrainRateT198R30.csv',1,0);
tension_strnT198R4     = csvread('tension_constTrueStrainRateT198R4.csv',1,0);
% Tension test under constant strain rate at temperature T = 298K
tension_strnT298R2000  = csvread('tension_constTrueStrainRateT298R2000.csv',1,0);
tension_strnT298R300   = csvread('tension_constTrueStrainRateT298R300.csv',1,0);
tension_strnT298R30    = csvread('tension_constTrueStrainRateT298R30.csv',1,0);
tension_strnT298R4     = csvread('tension_constTrueStrainRateT298R4.csv',1,0);
% Tension test under constant strain rate at temperature T = 398K
tension_strnT398R2000  = csvread('tension_constTrueStrainRateT398R2000.csv',1,0);
tension_strnT398R300   = csvread('tension_constTrueStrainRateT398R300.csv',1,0);
tension_strnT398R30    = csvread('tension_constTrueStrainRateT398R30.csv',1,0);
tension_strnT398R4     = csvread('tension_constTrueStrainRateT398R4.csv',1,0);

%% Data plot

% Plot true-stress vs. true-strain along tension direction at different 
% strain rate condition at temperatue T = 298K
% Strain in unit %, stress in unit MPa
ifg = ifg + 1; 
figure(ifg);
% Tension
plot(tension_strnT298R2000(:,11)*100,tension_strnT298R2000(:,12)/10^6,'-ob'); 
hold on
plot(tension_strnT298R300(:,11)*100,tension_strnT298R300(:,12)/10^6,'-xr'); 
hold on
plot(tension_strnT298R30(:,11)*100,tension_strnT298R30(:,12)/10^6,'-ok'); 
hold on
plot(tension_strnT298R4(:,11)*100,tension_strnT298R4(:,12)/10^6,'-.'); 
hold off
axis([-1,26,-0.1,2]);
title('Constant Strain rate at T=298K');
legend('2x10^{-2}','3x10^{-3}','3x10^{-4}','4x10^{-5}');
xlabel('\epsilon_{true} (%)');
ylabel('\sigma_{true} (MPa)');

% Plot true-stress vs. true-strain along tension direction under constant 
% strain rate condition epsilon_dt = 4x10^-5 at different temperatue
% Strain in unit %, stress in unit MPa
ifg = ifg + 1; 
figure(ifg);
% Tension
plot(tension_strnT198R4(:,11)*100,tension_strnT198R4(:,12)/10^6,'-ob'); 
hold on
plot(tension_strnT298R4(:,11)*100,tension_strnT298R4(:,12)/10^6,'-xr'); 
hold on
plot(tension_strnT398R4(:,11)*100,tension_strnT398R4(:,12)/10^6,'-ok'); 
hold off
axis([-1,26,-0.1,2]);
title('Constant Strain rate at d\epsilon/dt = 4x10^{-5}');
legend('T=198K','T=298K','T=398K');
xlabel('\epsilon_{true} (%)');
ylabel('\sigma_{true} (MPa)');