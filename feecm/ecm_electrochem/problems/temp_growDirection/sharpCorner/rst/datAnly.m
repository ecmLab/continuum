clc; clear; myFigureSetting;

%% parameter
% Geometry of the model
str = 'smth';
xSE   = 0.5;    % length of the SE, in unit um
ySE   = 0.5;    % height of the SE, in unit um
if strcmp(str,'mdlL1')
    xLi   = 0.1;   % length of the Li, in unit um
else
    xLi   = 0.02;   % length of the Li, in unit um
end
if strcmp(str,'mdlW1')
    yLi   = 0.01;   % height of the Li, in unit um
else
    yLi   = 0.001;  % height of the Li, in unit um
end
% Electrochemistry parameters
i_exc = 13;       % The exchange current density of lps, mA/cm^2
alfa  = 0.5;       % Reaction rate
i_applied = 1.0;  % The applied current density
F_RT  = 38.68; % Constant: F/RT, unit 1/V
ifg   = 0;      % Figure plot flag

%% Load data
% Total current at each boundary, in unit nA
% Column 1: Total current flow through Li_anode;   Column 2: Total current flow through SE_anode;
% Column 3: Total current flow through SE_cathode; Column 4: Total current flow through interface;
% positive means flow out from the block, negative means flow into the block.
Tot_crnt     = csvread([str,'.csv'],1,1);
% SE potential and current
SEAnode_pot  = csvread([str,'_SEAnode_potential_0001.csv'],1,0);
SEInt_pot    = csvread([str,'_SEInt_potential_0001.csv'],1,0);
SESide_pot   = csvread([str,'_SESide_potential_0001.csv'],1,0);
SEAnode_crnt = csvread([str,'_SEAnode_current_0001.csv'],1,0);
SEInt_crnt   = csvread([str,'_SEInt_current_0001.csv'],1,0);
% Li potential and current
LiCathode_pot  = csvread([str,'_LiCathode_potential_0001.csv'],1,0);
LiInt_pot    = csvread([str,'_LiInt_potential_0001.csv'],1,0);
LiSide_pot   = csvread([str,'_LiSide_potential_0001.csv'],1,0);
LiAnode_crnt = csvread([str,'_LiAnode_current_0001.csv'],1,0);

%% Data process
% Analytic solution of deposition rate according B-V equation
eta = SEInt_pot(:,5);  % Overpotential, the Li metal voltage is very close to zero (10^-10)
depRtInt = i_exc * (exp(alfa*F_RT*eta) - exp(-(1-alfa)*F_RT*eta));

% Percentage of electrons flow into the Li metal block to the total applied current
Tot_crnt(1)/Tot_crnt(3)

%% Plot
% Plot the potential of SE along thickness
ifg = ifg + 1;
figure(ifg)
plot(SEInt_pot(:,1),SEInt_pot(:,5),'.b',SESide_pot(:,1),SESide_pot(:,5),'-k')
legend('SE-pot-Top','SE-pot-Bottom');
% Plot the potential of SE at the anode
ifg = ifg + 1;
figure(ifg)
plot(SEAnode_pot(:,2),SEAnode_pot(:,5),'-k')
legend('SE-pot-Anode');
% Plot the normal current of SE at the anode
ifg = ifg + 1;
figure(ifg)
plot(SEAnode_crnt(:,2),-SEAnode_crnt(:,5),'-k')
legend('SE-crnt-Anode');
% Plot the normal current of SE at the interface, compare with analytic result
ifg = ifg + 1;
figure(ifg)
plot(SEInt_crnt(:,1),SEInt_crnt(:,6),'-k', SEInt_pot(:,1), depRtInt, '.b')
legend('SE-crnt-Interface','Analytic');

% Plot the potential of Li along thickness
ifg = ifg + 1;
figure(ifg)
plot(LiInt_pot(:,1),LiInt_pot(:,5),'o',LiSide_pot(:,1),LiSide_pot(:,5),'-k')
legend('Li-pot-Bottom','Li-pot-Top');
% Plot the potential of Li at the cathode
ifg = ifg + 1;
figure(ifg)
plot(LiCathode_pot(:,2),LiCathode_pot(:,5),'-k')
legend('Li-pot-Cathode');
% Plot the normal current of Li at the anode
ifg = ifg + 1;
figure(ifg)
plot(LiAnode_crnt(:,2),LiAnode_crnt(:,5),'-k')
legend('Li-crnt-Anode');
