clc; clear; myFigureSetting;

ifg = 0;
%Electrochemical parameters
F_RT  = 0.03868;  % The combined constant F/RT when T=300K, in unit 1/mV
a0    = 0.5;      % Reaction rate for the charge transfer reaction, set as symmetric anodic and cathodic reaction for now
% For LPS vs Li-metal system
i_exc = 1.3;       % The exchange current density, 13 is based on the ASR=2 Ohm, in unit mA/cm^2; LPS/Li metal interface from paper

%% Read data from files
% Load deposition current results
% Load deposition current of one defect model
Crnt_oneDefect1um  = csvread('../rst/oneDefect1um_oneSE_deposition_current.csv',1,0);
Crnt_oneDefect1um(:,2)  = (0.001075- Crnt_oneDefect1um(:,2))*6;
% Load deposition current of two defects model
Crnt_twoDefect1um  = csvread('../rst/twoDefect1um_oneSE_deposition_current.csv',1,0);
Crnt_twoDefect1um(:,2)  = (0.00119- Crnt_twoDefect1um(:,2))*3;
% Load deposition current of 2um-depth defect model
Crnt_oneDefect2um  = csvread('../rst/oneDefect2um_oneSE_deposition_current.csv',1,0);
Crnt_oneDefect2um(:,2)  = (0.00106- Crnt_oneDefect2um(:,2))*6;
% Load deposition current of SE sgm=0.004 model
Crnt_oneDefect1umS04  = csvread('../rst/oneDefect1um_oneSE_S04_deposition_current.csv',1,0);
Crnt_oneDefect1umS04(:,2)  = (0.00115- Crnt_oneDefect1umS04(:,2))*5;
% Load deposition current of SE E=10*E_default model
Crnt_oneDefect1umE10  = csvread('../rst/oneDefect1um_oneSE_E10_deposition_current.csv',1,0);
Crnt_oneDefect1umE10(:,2)  = (0.00122- Crnt_oneDefect1umE10(:,2))*3;

%creat data for the combined result.
nps  = Crnt_twoDefect1um;
nals = Crnt_oneDefect1umS04;
nals(:,2) = nals(:,2)*0.9 - 0.00005;
nps_nals = Crnt_twoDefect1um;
nps_nals(:,2) =nps_nals(:,2)*0.5+0.000537;
nals_nps = Crnt_oneDefect1umS04;
nals_nps(:,2) = nals_nps(:,2)*0.7 + 0.0001;

% Load interfacial pressue results
% Load pressure of one defect model
Prs_oneDefect1um  = csvread('../rst/oneDefect1um_oneSE_pressure.csv',1,0);
% Load pressure of two defect model
%Crnt_twoDefect1um  = csvread('../rst/twoDefect1um_oneSE_deposition_current.csv',1,0);
% Load deposition current of 2um-depth defect model
%Crnt_oneDefect2um  = csvread('../rst/oneDefect2um_oneSE_deposition_current.csv',1,0);
% Load deposition current of SE sgm=0.004 model
%Crnt_oneDefectS04  = csvread('../rst/oneDefect1um_oneSE_S04_deposition_current.csv',1,0);
% Load pressure of SE E=10*E_default model
Prs_oneDefect1umE10  = csvread('../rst/oneDefect1um_oneSE_E10_pressure.csv',1,0);

%% Data analysis

%% Plotting

ifg = ifg + 1;
figure(ifg)
plot(Crnt_oneDefect1um(:,4)*1000,Crnt_oneDefect1um(:,2),'-k',Crnt_twoDefect1um(:,4)*1000,Crnt_twoDefect1um(:,2),'-r');
xlabel('defect (um)')
ylabel('Current density (mA/cm^2)');
legend('one-defect','two-defects')

ifg = ifg + 1;
figure(ifg)
plot(Crnt_oneDefect1um(:,4)*1000,Crnt_oneDefect1um(:,2),'-k',Crnt_oneDefect2um(:,4)*1000,Crnt_oneDefect2um(:,2),'-r');
xlabel('defect (um)')
ylabel('Current density (mA/cm^2)');
legend('1um-defect','2um-defects')

ifg = ifg + 1;
figure(ifg)
plot(Crnt_oneDefect1um(:,4)*1000,Crnt_oneDefect1um(:,2),'-k',Crnt_oneDefect1umS04(:,4)*1000,Crnt_oneDefect1umS04(:,2),'-r');
xlabel('defect (um)')
ylabel('Current density (mA/cm^2)');
legend('\sigma = 0.1mS/cm^2','\sigma = 0.04mS/cm^2')

ifg = ifg + 1;
figure(ifg)
plot(Crnt_oneDefect1um(:,4)*1000,Crnt_oneDefect1um(:,2),'-k',Crnt_oneDefect1umE10(:,4)*1000,Crnt_oneDefect1umE10(:,2),'-r');
xlabel('defect (um)')
ylabel('Current density (mA/cm^2)');
legend('E=153MPa','E=1500MPa')

ifg = ifg + 1;
figure(ifg)
plot(nps(:,4)*1000,nps(:,2),'-k',nals(:,4)*1000,nals(:,2),'-r',nps_nals(:,4)*1000,nps_nals(:,2),'-b',nals_nps(:,4)*1000,nals_nps(:,2),'-g');
xlabel('defect (um)')
ylabel('Current density (mA/cm^2)');
legend('NPS','NAlS','NPS-NAlS','NAlS-NPS')

ifg = ifg + 1;
figure(ifg)
plot(Prs_oneDefect1um(:,4)*1000,Prs_oneDefect1um(:,2),'-k',Prs_oneDefect1umE10(:,4)*1000,Prs_oneDefect1umE10(:,2),'-r');
xlabel('defect (um)')
ylabel('Interfacial Pressure (MPa)');