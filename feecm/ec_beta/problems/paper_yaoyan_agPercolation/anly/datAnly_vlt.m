clc; clear; myFigureSetting;

ifg = 0;
%Geometric parameter
xBL = 10;      % width of the BL, in unit um
yBL = 10;      % thickness of the BL, in unit um
dia = 1.0;      % diameter of pore and Ag, in unit um
%Electrochemical parameters
F_RT  = 0.03868;  % The combined constant F/RT when T=300K, in unit 1/mV
a0    = 0.5;      % Reaction rate for the charge transfer reaction, set as symmetric anodic and cathodic reaction for now
% For LPS vs Li-metal system
iSE_exc = 1.3;        % The exchange current density, 13 is based on the ASR=2 Ohm, in unit mA/cm^2; LPS/Li metal interface from paper
iAg_exc = iSE_exc*5;  % The exchange current density, 13 is based on the ASR=2 Ohm, in unit mA/cm^2; LPS/Li metal interface from paper
iCC_exc = iSE_exc*10; % The exchange current density, 13 is based on the ASR=2 Ohm, in unit mA/cm^2; LPS/Li metal interface from paper
% For change variables
sgmLi = 10.^(linspace(-1.1,1.1,12));  % The ionic conductivity, set as variable for parameter studies, in unit mS/cm
sgmEn = 10.^(linspace(1.0,3.6,14));  % The ionic conductivity, set as variable for parameter studies, in unit mS/cm
vAg   = linspace(-1,-10,10); % The equilibrium voltage of the AgLi particle in the surface defect, in unit mV; 0 corresponds to Li metal anode;

for jAg = 1 : length(vAg)
%% Read data from files
% Relation of filenames and variable names: "sgm001" and "sgm-2" means sgm = 10^-2 = 0.01 mS/cm 
% Load potentials at different interfaces
    tmpStr  = ['sePot_Ag',num2str(jAg),'= csvread(''../rst/Vlt_sePot_Ag-',num2str(jAg),'.csv'',1,0);'];
    eval(tmpStr);
    tmpStr  = ['porePot_Ag',num2str(jAg),'= csvread(''../rst/Vlt_porePot_Ag-',num2str(jAg),'.csv'',1,0);'];
    eval(tmpStr);
    tmpStr  = ['agPot_Ag',num2str(jAg),'= csvread(''../rst/Vlt_agPot_Ag-',num2str(jAg),'.csv'',1,0);'];
    eval(tmpStr);
    tmpStr  = ['ccPot_Ag',num2str(jAg),'= csvread(''../rst/Vlt_ccPot_Ag-',num2str(jAg),'.csv'',1,0);'];
    eval(tmpStr);
    tmpStr  = ['x1Pot_Ag',num2str(jAg),'= csvread(''../rst/Vlt_z1Pot_Ag-',num2str(jAg),'.csv'',1,0);'];
    eval(tmpStr);
    tmpStr  = ['x2Pot_Ag',num2str(jAg),'= csvread(''../rst/Vlt_z2Pot_Ag-',num2str(jAg),'.csv'',1,0);'];
    eval(tmpStr);
    tmpStr  = ['x3Pot_Ag',num2str(jAg),'= csvread(''../rst/Vlt_z3Pot_Ag-',num2str(jAg),'.csv'',1,0);'];
    eval(tmpStr);

% porePot_Ag1  = csvread('../rst/mdl_pore_potential_0001.csv',1,0);
% agPot_Ag1    = csvread('../rst/mdl_ag_potential_0001.csv',1,0);
% ccPot_Ag1    = csvread('../rst/mdl_CC_potential_0001.csv',1,0);
% x1Pot_Ag1    = csvread('../rst/mdl_zeroX1_potential_0001.csv',1,0);
% x2Pot_Ag1    = csvread('../rst/mdl_zeroX2_potential_0001.csv',1,0);
% x3Pot_Ag1    = csvread('../rst/mdl_zeroX3_potential_0001.csv',1,0);

%% Data analysis
% Compute the Li potential along the middle line
    tmpStr  = ['midPot_Ag',num2str(jAg),'= [x1Pot_Ag',num2str(jAg),';porePot_Ag',num2str(jAg),'(porePot_Ag',...
                num2str(jAg),'(:,5)>0,:);x2Pot_Ag',num2str(jAg),';agPot_Ag',num2str(jAg),'(agPot_Ag',num2str(jAg),'(:,5)>0,:);x3Pot_Ag',num2str(jAg),'];'];
    eval(tmpStr);
% Compute normal current density from the analytic Butler-Volmer relation:
 % i_n = i_exc * (exp(alpha*eta*F/RT)-exp(-alpha*eta*F/RT))
    tmpStr  = ['seEta_Ag',num2str(jAg),'  = - sePot_Ag',num2str(jAg),'(:,2) - sePot_Ag',num2str(jAg),'(:,3);'];
    eval(tmpStr);
    tmpStr  = ['seCrnt_Ag',num2str(jAg),' = iSE_exc * (exp(a0*F_RT*seEta_Ag',num2str(jAg),') - exp(-a0*F_RT*seEta_Ag',num2str(jAg),'));'];
    eval(tmpStr);
    tmpStr  = ['poreEta_Ag',num2str(jAg),'  = - porePot_Ag',num2str(jAg),'(:,2) - porePot_Ag',num2str(jAg),'(:,3);']; 
    eval(tmpStr);
    tmpStr  = ['poreCrnt_Ag',num2str(jAg),' = iSE_exc * (exp(a0*F_RT*poreEta_Ag',num2str(jAg),') - exp(-a0*F_RT*poreEta_Ag',num2str(jAg),'));'];
    eval(tmpStr);
    tmpStr  = ['agEta_Ag',num2str(jAg),'  = vAg(',num2str(jAg),') - agPot_Ag',num2str(jAg),'(:,2) - agPot_Ag',num2str(jAg),'(:,3);']; 
    eval(tmpStr);
    tmpStr  = ['agCrnt_Ag',num2str(jAg),' = iAg_exc * (exp(a0*F_RT*agEta_Ag',num2str(jAg),') - exp(-a0*F_RT*agEta_Ag',num2str(jAg),'));'];
    eval(tmpStr);
    tmpStr  = ['ccEta_Ag',num2str(jAg),'  = - ccPot_Ag',num2str(jAg),'(:,2) - ccPot_Ag',num2str(jAg),'(:,3);'];
    eval(tmpStr);
    tmpStr  = ['ccCrnt_Ag',num2str(jAg),' = iCC_exc * (exp(a0*F_RT*ccEta_Ag',num2str(jAg),') - exp(-a0*F_RT*ccEta_Ag',num2str(jAg),'));'];
    eval(tmpStr);

% Compute the total current ratio at each interfaces
    tmpStr   =['seIave(jAg) = sum(seCrnt_Ag',num2str(jAg),')/length(seCrnt_Ag',num2str(jAg),');'];  eval(tmpStr);
    tmpStr   =['ccIave(jAg) = sum(ccCrnt_Ag',num2str(jAg),')/length(ccCrnt_Ag',num2str(jAg),');'];  eval(tmpStr);
    tmpStr   =['poreIave(jAg) = sum(poreCrnt_Ag',num2str(jAg),')/length(poreCrnt_Ag',num2str(jAg),');'];  eval(tmpStr);
    tmpStr   =['agIave(jAg) = sum(agCrnt_Ag',num2str(jAg),')/length(agCrnt_Ag',num2str(jAg),');'];  eval(tmpStr);

    % Compute the total current ratio at each interfaces
    tmpStr   =['seIr(jAg) = -mean(seCrnt_Ag',num2str(jAg),');'];  eval(tmpStr);
    tmpStr   =['ccIr(jAg) = -mean(ccCrnt_Ag',num2str(jAg),');'];  eval(tmpStr);
    tmpStr   =['poreIr(jAg) = -mean(poreCrnt_Ag',num2str(jAg),')*pi*dia/yBL;'];  eval(tmpStr);
    tmpStr   =['agIr(jAg) = -mean(agCrnt_Ag',num2str(jAg),')*pi*dia/yBL;'];  eval(tmpStr);

% Compute the Li potential along the middle line
% midPot_Ag1   = [x1Pot_Ag1;porePot_Ag1(porePot_Ag1(:,5)>0,:);x2Pot_Ag1;agPot_Ag1(agPot_Ag1(:,5)>0,:);x3Pot_Ag1];

% % Compute normal current density from the analytic Butler-Volmer relation:
%  % i_n = i_exc * (exp(alpha*eta*F/RT)-exp(-alpha*eta*F/RT))
% seEta_Ag1  = - sePot_Ag1(:,2) - sePot_Ag1(:,3); 
% seCrnt_Ag1 = iSE_exc * (exp(a0*F_RT*seEta_Ag1) - exp(-a0*F_RT*seEta_Ag1));
% poreEta_Ag1  = - porePot_Ag1(:,2) - porePot_Ag1(:,3); 
% poreCrnt_Ag1 = iSE_exc * (exp(a0*F_RT*poreEta_Ag1) - exp(-a0*F_RT*poreEta_Ag1));
% agEta_Ag1  = vAg(1) - agPot_Ag1(:,2) - agPot_Ag1(:,3); 
% agCrnt_Ag1 = iAg_exc * (exp(a0*F_RT*agEta_Ag1) - exp(-a0*F_RT*agEta_Ag1));
% ccEta_Ag1  = - ccPot_Ag1(:,2) - ccPot_Ag1(:,3); 
% ccCrnt_Ag1 = iCC_exc * (exp(a0*F_RT*ccEta_Ag1) - exp(-a0*F_RT*ccEta_Ag1));

end

%% Plotting
% Plot potential at sgm= 0.8 mS/cm, ASR= 2 Ohm*cm^2 as a function of x coordinate at four different interfaces
ifg = ifg + 1;
figure(ifg)
plot(sePot_Ag1(:,4),sePot_Ag1(:,3),'-b',ccPot_Ag1(:,4),ccPot_Ag1(:,3),'-r',midPot_Ag1(:,4),midPot_Ag1(:,3),'-k');
legend('SE','CC','Center');
title('Voltage at Interface, in mV');
% Plot the current at sgm= 0.1 mS/cm, ASR= 2 Ohm*cm^2 and varying AgLi particle equilibrium potential in the surface defect = [0, -1, -2, -3, -4 ,-5] mV
ifg = ifg + 1;
figure(ifg)
plot(sePot_Ag1(:,4),-seCrnt_Ag1,'-b',ccPot_Ag1(:,4),-ccCrnt_Ag1,'-r',...
     porePot_Ag1(porePot_Ag1(:,5)>0,4),-poreCrnt_Ag1(porePot_Ag1(:,5)>0),'-k',agPot_Ag1(porePot_Ag1(:,5)>0,4),-agCrnt_Ag1(porePot_Ag1(:,5)>0),'-k');
legend('SE','CC','Middle');
title('Current at SE/Anode Interface, in mA/cm^2');

% Plot the current at sgm= 0.1 mS/cm, ASR= 2 Ohm*cm^2 and varying AgLi particle equilibrium potential in the surface defect = [0, -1, -2, -3, -4 ,-5] mV
ifg = ifg + 1;
figure(ifg)
agLiy = [5.212406949875933115, 9.400812534416713717, 10.92401456535275805, 12.4336165781554282, 14.13881885175847053,...
         16.21202161602882796, 18.91562522083362552, 22.77083036110713365, 29.11083881445175692,43.21165761554355811];
plot(agLiy,flip(seIr)*100,'-b',agLiy,flip(ccIr)*100,'-r',agLiy,flip(poreIr)*100,'-k',agLiy,flip(agIr)*100,'-y');
% plot(-vAg,seIr*100,'-b',-vAg,ccIr*100,'-r',-vAg,poreIr*100,'-k',-vAg,agIr*100,'-y');
legend('SE','CC','Pore','Ag');
title('Current ratio (%)');
