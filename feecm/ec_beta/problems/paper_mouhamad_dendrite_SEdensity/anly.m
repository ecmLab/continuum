clc; clear;
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontSize',16);

%% parameter
% Universal constants
F_RT  = 38.68;         % Constant: F/RT, unit 1/V
VLi_F = 1.347*10^-4;   % MolarVolume_Li/F, in unit cm^3/(s*A)
ifg   = 0;             % Figure plot flag
% Electrochemistry parameters
sgmSE = 0.1;       % The ionic conductivity of the SE, mS/cm
% Geometric parameters
hSE   = 20.0;  % thickness of the SE, in unit um
wSE   = 8.0;   % width of the SE, in unit um
rp    = 0.05;  % the round radius at corner point, in unit um
% Temporal parameters
dt    = 0.1;       % Timestep, in unit hours
% Load data
load hLi.mat;
load tLi.mat;
% Variables
nT    = size(tLi,2);       % Total number of steps
nSmp  = size(tLi,1);        % Number of sampling points for recoording Li thickness
yLi   = linspace(0,hSE-1.5*rp,nSmp)'; % y Coordinate of sampling points, ending at 1.5*rp below top surface

%% Data analysis
ifg = ifg + 1;
figure(ifg)
for i = 1:nT
    tplt = tLi(tLi(:,i)>0,i);
    yplt = yLi(tLi(:,i)>0);
    plot(yplt,tplt,'.')
    hold on
end
hold off

ifg = ifg + 1;
figure(ifg)
for it = 1 : nT
% Load CSV data from MOOSE calculation: structure of the csv data
% Column 1-2: current density along the x and y direction, in unit mA/cm^2
% Column 4-5: coordinate in the x and y direction, in unit um
    tmp         = csvread(['rst/t',num2str(it),'.csv'],1,0);
    [cord,indx] = sortrows(tmp(:,4:5),[1,2],{'ascend' 'descend'}); % sort the coordinate with x ascend and y dscend
    crnt        = tmp(indx,1:2); % Sort the current density accordingly
% Compute the arch length along interface, and the normal current density    
    archL    = zeros(length(cord),1);
    nrCnt    = zeros(length(cord),1);
    nrCnt(1) = crnt(1,2);
    for ip = 2 : length(cord)
        archL(ip) = archL(ip-1) + norm(cord(ip,:)-cord(ip-1,:));
        if cord(ip,1) == wSE
            nrCnt(ip) = crnt(ip,1);
        elseif cord(ip,2) == hSE
            nrCnt(ip) = crnt(ip,2);
        else
            dNr       = cord(ip,:) - [wSE-rp,hSE-rp];
            nrCnt(ip) = crnt(ip,:)*dNr'/norm(dNr);
        end
    end
% Plot
    plot(archL, nrCnt)
    hold on
end
hold off
% axis([15,20,0,1.5])
%% Plot

