%% This file include all the equilibrium voltage of different material as a function of Li content
clc;clear;
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontName','Arial');
set(0,'defaultAxesFontSize',18);
ifg         = 0;

%% Read data
FCC       = csvread('FCC.csv',1,0);       % 1st column is x in Ag_(1-x)Li_x; 2nd column is free energy, in unit kJ/mol
beta      = csvread('beta.csv',1,0);      % 1st column is x in Ag_(1-x)Li_x; 2nd column is free energy, in unit kJ/mol
gama      = csvread('gamma.csv',1,0);     % 1st column is x in Ag_(1-x)Li_x; 2nd column is free energy, in unit kJ/mol
BCC       = csvread('BCC.csv',1,0);       % 1st column is x in Ag_(1-x)Li_x; 2nd column is free energy, in unit kJ/mol

%% Data analysis
% 1. Fit with fit function
FCCq      = fit(FCC(:,1),FCC(:,2),'smoothingspline');
betaq     = fit(beta(:,1),beta(:,2),'poly5');
gamaq     = fit(gama(:,1),gama(:,2),'poly2');
BCCq      = fit(BCC(:,1),BCC(:,2),'poly1');

% 2. plot with Fitted data
% Extract parameter from the fitted function, assuming a five order polynomial: y = p1 + p2*x + p3*x^2 + p4*x^3 + p5*x^4 + p6*x^5
% 1th row: parameters for FCC; 2nd row: parameters for beta; 3rd row:parameters for gama; 4th row: parameters for bcc
AA    = [0.1825, -61.39, 59.00, 0.0, 0.0, 0.0; 1916.9, -14062,  40108, -56284, 38879, -10544;...
         29.53,  -164.4, 141.1, 0.0, 0.0, 0.0; -51.169, 50.806, 0.0,   0.0,    0.0,   0.0];
x     = linspace(0.0,1,100);
FCf         =  AA(1,1) + AA(1,2)*x + AA(1,3)*x.^2 + AA(1,4)*x.^3 + AA(1,5)*x.^4 + AA(1,6)*x.^5;
FCf(x>0.75)  = NaN;
BTf         =  AA(2,1) + AA(2,2)*x + AA(2,3)*x.^2 + AA(2,4)*x.^3 + AA(2,5)*x.^4 + AA(2,6)*x.^5;
BTf(x<0.33) = NaN;
GMf         =  AA(3,1) + AA(3,2)*x + AA(3,3)*x.^2 + AA(3,4)*x.^3 + AA(3,5)*x.^4 + AA(3,6)*x.^5;
BCf         =  AA(4,1) + AA(4,2)*x + AA(4,3)*x.^2 + AA(4,4)*x.^3 + AA(4,5)*x.^4 + AA(4,6)*x.^5;
BCf(x<0.75)  = NaN;

% 3. Find common tangents points between two phases
[FB_P1, FB_P2]   = comTangent(AA(1,:),AA(2,:),[0.01;0.3]);   % Common tangents between FCC and beta
[BG_P1, BG_P2]   = comTangent(AA(2,:),AA(3,:),[0.5;0.9]);    % Common tangents between beta and gamma

%% Plot
% 1. Plot voltage as a function of y in MLi_y
ifg = ifg + 1;
figure(ifg)
plot(FCC(:,1),FCC(:,2),'.r', beta(:,1),beta(:,2),'.r', gama(:,1),gama(:,2),'.r', BCC(:,1),BCC(:,2),'.r')
hold on
plot(x,FCf,'b', x,BTf,'b', x,GMf,'b', x,BCf,'b')
hold on
plot([FB_P1(1,1),FB_P2(1,1)],[FB_P1(1,2),FB_P2(1,2)],'--ok')
hold on
plot([BG_P1(1,1),BG_P2(1,1)],[BG_P1(1,2),BG_P2(1,2)],'--ok',[BG_P1(2,1),BG_P2(2,1)],[BG_P1(2,2),BG_P2(2,2)],'--ok')
hold off