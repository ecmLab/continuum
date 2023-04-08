%% This file include all the nucleation overpotential of different metal material
clc;clear;

%% Parameters
i_t       = 0.01;               % applied current density, unit mA/cm^2

%% Read data
% 1st column is capacity; 2nd column is voltage vs Li/Li+, in unit V
AuLi      = csvread('nucleation_Au.csv',0,0);
AgLi      = csvread('nucleation_Ag.csv',0,0);
ZnLi      = csvread('nucleation_Zn.csv',0,0);
MgLi      = csvread('nucleation_Mg.csv',0,0);
AlLi      = csvread('nucleation_Al.csv',0,0);
PtLi      = csvread('nucleation_Pt.csv',0,0);
CuLi      = csvread('nucleation_Cu.csv',0,0);
NiLi      = csvread('nucleation_Ni.csv',0,0);
CLi       = csvread('nucleation_C.csv',0,0);
SnLi      = csvread('nucleation_Sn.csv',0,0);
SiLi      = csvread('nucleation_Si.csv',0,0);

%% Data analysis
AuLi(:,2) = AuLi(:,2) - 0.05;    % Voltage of Au is shift by 0.05V
AgLi(:,2) = AgLi(:,2) - 0.10;    % Voltage of Ag is shift by 0.10V
ZnLi(:,2) = ZnLi(:,2) - 0.15;    % Voltage of Zn is shift by 0.15V
MgLi(:,2) = MgLi(:,2) - 0.20;    % Voltage of Mg is shift by 0.20V
AlLi(:,2) = AlLi(:,2) - 0.25;    % Voltage of Al is shift by 0.25V
PtLi(:,2) = PtLi(:,2) - 0.30;    % Voltage of Pt is shift by 0.30V
CuLi(:,2) = CuLi(:,2) - 0.05;    % Voltage of Cu is shift by 0.05V
NiLi(:,2) = NiLi(:,2) - 0.10;    % Voltage of Ni is shift by 0.10V
CLi(:,2)  = CLi(:,2)  - 0.15;    % Voltage of C  is shift by 0.15V
SnLi(:,2) = SnLi(:,2) - 0.20;    % Voltage of Sn is shift by 0.20V
SiLi(:,2) = SiLi(:,2) - 0.25;    % Voltage of Si is shift by 0.25V

depV   = [AuLi(end,2),AgLi(end,2),ZnLi(end,2),MgLi(end,2),AlLi(end,2),PtLi(end,2),...
          CuLi(end,2),NiLi(end,2),CLi(end,2), SnLi(end,2),SiLi(end,2)]'*1000;  % The overpotential of Li plating, in unit mV
R_ct   = depV/i_t;        % ASR, in unit Ohm*cm^2
N_en   = [min(AuLi(:,2)),min(AgLi(:,2)),min(ZnLi(:,2)),min(MgLi(:,2)),min(AlLi(:,2)),min(PtLi(:,2)),...
          min(CuLi(:,2)),min(NiLi(:,2)),min(CLi(:,2)), min(SnLi(:,2)),min(SiLi(:,2))]'*1000 - i_t*R_ct; % nucleation overpotential, in unit mV

%% Plot
ifg = 0;
% 1. Plot voltage as a function of capacity
ifg = ifg + 1;
figure(ifg)
plot(AuLi(:,1),AuLi(:,2),'-r', AgLi(:,1),AgLi(:,2),'-y', ZnLi(:,1),ZnLi(:,2),'-g', SnLi(:,1),SnLi(:,2),'-p',...
     MgLi(:,1),MgLi(:,2),'-p', AlLi(:,1),AlLi(:,2),'-p', PtLi(:,1),PtLi(:,2),'-p',CLi(:,1),CLi(:,2),'-g');
legend('Au','Ag','Zn','Sn', 'Mg','Al', 'Pt', 'C');
ifg = ifg + 1;
figure(ifg)
plot(CuLi(:,1),CuLi(:,2),'-r', NiLi(:,1),NiLi(:,2),'-y', CLi(:,1),CLi(:,2),'-g', ...
     SnLi(:,1),SnLi(:,2),'-p', SiLi(:,1),SiLi(:,2),'-p');
legend('Cu','Ni','C','Sn','Si');
 
 % 2. Plot ASR, in unit Ohm*cm^2
ifg = ifg + 1;
figure(ifg)
plot(1:length(R_ct), -R_ct, 'o')

% 2. Plot nucleation energy, in unit mV
ifg = ifg + 1;
figure(ifg)
plot(1:length(R_ct), -N_en, 'o')