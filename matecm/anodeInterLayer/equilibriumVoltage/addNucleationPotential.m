%% This file include the equilibrium voltage and nucleation overpotential of Sn
clc;clear;
FF          = 96485.33;              % Faraday constant, in unit s*A/mol

%% Read data
% Equilibrium voltage of the Li-Sn system
SnLi        = csvread('SnLi.csv',1,0);    % 1st column is x in SnLix; 2nd column is voltage vs Li/Li+, in unit V
SnLi        = SnLi(1:end-1,:);
SnLi(end,1) = 4.4;
SnLi(:,2)   = SnLi(:,2)*1000;              % Convert the 2nd column to the unit of mV
% Nucleation voltage of the Li-Sn system
SnLinuc      = csvread('nucleation_Sn.csv',0,0);
SnLinuc(:,2) = SnLinuc(:,2) - 0.20;    % Voltage of Sn is shift by 0.20V

%% Data analysis
% 1. Convert x coordinate of nucleation from capacity to Li content
% 2. Add the x coordinate of nucleation to the end of equilibrium voltage, with 2~5 shift
SnLinuc(:,1) = SnLinuc(:,1)*100 + SnLi(end,1);
SnLinuc(:,2) = (SnLinuc(:,2) - SnLinuc(end,2))*1000;    % Voltage normalize by the end, and convert to mV

totV   =  [[0,SnLi(1,2)];SnLi;SnLinuc];
tV     =  smooth(totV(:,1),totV(:,2));
tV     = tV - tV(end);

%% Save
% 1. Save the equilibrium voltage for future usage
chemP_SnLi_db=[totV(:,1),-tV*FF/1000]; % Convert to J/mol
save('chemP_SnLi_db.mat','chemP_SnLi_db')

%% Plot
ifg = 0;
% 1. Plot voltage as a function of y in MLi_y
ifg = ifg + 1;
figure(ifg)
plot(SnLi(:,1),SnLi(:,2),'b', SnLinuc(:,1),SnLinuc(:,2))
ifg = ifg + 1;
figure(ifg)
plot(totV(:,1),totV(:,2),'b',totV(:,1),tV)