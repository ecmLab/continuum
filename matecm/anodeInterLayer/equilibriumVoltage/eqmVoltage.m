%% This file include all the equilibrium voltage of different material as a function of Li content
clc;clear;

%% Read data
% For Li-C system
CLi       = csvread('CLi.csv',1,0);     % 1st column is x in LiC(6/x); 2nd column is voltage vs Li/Li+, in unit V
CLi(:,1)  = CLi(:,1)/6;                 % Convert the 1st column to the format of CLiy
CLi(:,2)  = CLi(:,2)*1000;              % Convert the 2nd column to the unit of mV
% For Li-Ag system
AgLi      = csvread('AgLi.csv',1,0);    % 1st column is x in AgLix; 2nd column is voltage vs Li/Li+, in unit mV
% For Li-Sn system
SnLi      = csvread('SnLi.csv',1,0);    % 1st column is x in SnLix; 2nd column is voltage vs Li/Li+, in unit V
SnLi(:,2) = SnLi(:,2)*1000;              % Convert the 2nd column to the unit of mV
% For Li-Al system
AlLi      = csvread('AlLi.csv',1,0);    % 1st column is x in AlLix; 2nd column is voltage vs Li/Li+, in unit V
AlLi(:,2) = AlLi(:,2)*1000;              % Convert the 2nd column to the unit of mV

%% Plot
ifg = 0;
% 1. Plot voltage as a function of y in MLi_y
ifg = ifg + 1;
figure(ifg)
plot(CLi(:,1),CLi(:,2),'--b', AgLi(:,1),AgLi(:,2),'r', SnLi(:,1),SnLi(:,2),'b', AlLi(:,1),AlLi(:,2),'g')

% 2. Plot voltage as a function of x in M_(1-x)Li_x
ifg = ifg + 1;
figure(ifg)
xCLi  = CLi(:,1)./(1+CLi(:,1));
xAgLi = AgLi(:,1)./(1+AgLi(:,1));
xSnLi = SnLi(:,1)./(1+SnLi(:,1));
xAlLi = AlLi(:,1)./(1+AlLi(:,1));
plot(xCLi,CLi(:,2),'--b', xAgLi,AgLi(:,2),'r', xSnLi,SnLi(:,2),'b', xAlLi,AlLi(:,2),'g')



SnLinuc      = csvread('nucleation_Sn.csv',0,0);
SnLinuc(:,2) = SnLinuc(:,2) - 0.20;    % Voltage of Sn is shift by 0.20V
plot(SnLi(:,1),SnLi(:,2),'b', SnLinuc(:,1)*1000,SnLinuc(:,2)*1000)