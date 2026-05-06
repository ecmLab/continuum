clc; clear;
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontSize',16);

%% parameter
% Universal constants
F_RT  = 0.03868;  % The combined constant F/RT when T=300K, in unit 1/mV
VLi_F = 1.347*10^-4;   % MolarVolume_Li/F, in unit cm^3/(s*A)
ifg   = 0;             % Figure plot flag
% Electrochemistry parameters
i_exc = 1.3;       % The exchange current density, 13 is based on the ASR=2 Ohm, in unit mA/cm^2;
sgmSE = 0.1;           % The ionic conductivity of the SE, mS/cm
% Geometric parameters
wSE   = 0.5;  % width of the SE, in unit um
hSE   = 10;  % thickness of the SE, in unit um
dw    = 0.05;      % width of the crack, in unit um
hLi   = 2;      % length of Li, in unit um
hAg   = 1;  % center location of the Ag, in unit um
dAg   = 0.05;    % radius of the Ag particle, in unit um
potLiAg = -2;   % The potential of LiAg, in unit mV

%% Load CSV data from MOOSE calculation: structure of the csv data
% load potential
% Column 2: potential along the li metal, in unit mA/cm^2
% Column 3-4: coordinate in the x and y direction, in unit um
Ag1Ref      = csvread(['rst/1AgRef_anode_potential_0001.csv'],1,0);
Ag1Ref      = sortrows(Ag1Ref,[4,3],{'ascend' 'descend'});       % sort the coordinate with Y and x ascend
cord        = Ag1Ref(:,3:4); % coordinates of points
% Compute the arch length along Li metal, and the normal current density
archL    = zeros(length(cord),1);
for ip = 2 : length(cord)
    archL(ip) = archL(ip-1) + norm(cord(ip,:)-cord(ip-1,:));
end
Ag1Ref = [archL, i_exc * (exp(0.5*F_RT*Ag1Ref(:,2)) - exp(-0.5*F_RT*Ag1Ref(:,2)))];

% load Ag-LPS results
Ag1      = csvread(['rst/1Ag_anode_potential_0001.csv'],1,0);
Ag1      = sortrows(Ag1,[4,3],{'ascend' 'descend'});       % sort the coordinate with Y and x ascend
cord     = Ag1(:,3:4); % coordinates of points
% Compute the arch length along Li metal, and the normal current density
archL    = zeros(length(cord),1);
for ip = 2 : length(cord)
    archL(ip) = archL(ip-1) + norm(cord(ip,:)-cord(ip-1,:));
%     if cord(ip,1)==dw && abs(cord(ip,2) - hAg)<=dAg/2
    if cord(ip,1)>dw && cord(ip,2)>0
       Ag1(ip,2) = Ag1(ip,2) - potLiAg;
    end
end
Ag1 = [archL, i_exc * (exp(0.5*F_RT*Ag1(:,2)) - exp(-0.5*F_RT*Ag1(:,2)))];

%% Data analysis
ifg = ifg + 1;
figure(ifg)
plot(Ag1Ref(:,1), abs(Ag1Ref(:,2)),'--k',Ag1(:,1), abs(Ag1(:,2)),'-b')
% plot(Ag1Ref(:,1), Ag1Ref(:,2),'--k')


