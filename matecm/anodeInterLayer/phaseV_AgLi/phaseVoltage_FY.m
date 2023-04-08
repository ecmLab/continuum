%% This file include all the equilibrium voltage of different material as a function of Li content
clc;clear;
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontName','Arial');
set(0,'defaultAxesFontSize',18);
ifg         = 0;

%% 1. Read spline data
% 1st row: x for BCC, 2nd row: y for BCC;
% 2st row: x for FCC, 2nd row: y for FCC;
% 3st row: x for gama, 2nd row: y for gama;
bsp     = csvread('voltage_bspline.csv',0,0);       
BCC     = bsp(1:2,:)';
FCC     = bsp(3:4,:)';
gma     = bsp(5:6,:)';

%% 2. read potential data
eqmV     = load('FY_equilibrium.txt');
ovBV     = load('FY_overshoot_bcc.txt');
ovFV     = load('FY_overshoot_fcc.txt');
ovGV     = load('FY_overshoot_gam.txt');

%% Plot
% 1. Plot voltage as a function of y in MLi_y
ifg = ifg + 1;
figure(ifg)
plot(eqmV(:,1),eqmV(:,2), ovBV(:,1),ovBV(:,2),ovFV(:,1),ovFV(:,2),ovGV(:,1),ovGV(:,2))
% 2. Plot voltage as a function of x in M_(1-x)Li_x
ifg = ifg + 1;
figure(ifg)
plot(eqmV(:,1)./(1+eqmV(:,1)),eqmV(:,2), ovBV(:,1)./(1+ovBV(:,1)),ovBV(:,2),...
     ovFV(:,1)./(1+ovFV(:,1)),ovFV(:,2),ovGV(:,1)./(1+ovGV(:,1)),ovGV(:,2));
legend('eqm','BCC-over','FCC-over','gama-over');