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
i_exc = 0.001;  % The exchange current density, mA/cm^2
alfa  = 0.5;  % Reaction rate
c_ref = 0.0528;

% Geometry of the model
lx    =  1e-3;                      % length of the model, in unit um
lyz   =  1e-3;                      % height of the model, in unit um

%% Load data
% Meaning of each column in the csv file:
% 1st: time, in unit second;                                            
% 2st: equilibrium vlotage;
% 3st: cell voltage;                     
% 4st: chemical potential
% 5st: concentration;
rst  = csvread('test_ag.csv',1,0);

dlt  = rst(:,3) -rst(:,2);
crnt = i_exc*FF*dlt/RT/1000;

%% Data plot
% Plot true-stress vs. time along tension direction 
% stress in unit MPa, time in unit seconds
ifg = ifg + 1; 
figure(ifg);
plot(rst(:,1)/3600,crnt); 

ifg = ifg + 1; 
figure(ifg);
plot(rst(:,1)/3600,dlt);

ifg = ifg + 1; 
figure(ifg);
plot(rst(:,5)/c_ref,rst(:,2),'-b',rst(:,5)/c_ref,rst(:,3),'-k');

