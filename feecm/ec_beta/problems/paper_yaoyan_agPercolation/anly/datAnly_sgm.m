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
sgmLi = 10.^[-1.1:0.2:2.0];  % The ionic conductivity, set as variable for parameter studies, in unit mS/cm
sgmEn = 10.^[0.01:0.2:3.01];  % The ionic conductivity, set as variable for parameter studies, in unit mS/cm
vAg   = linspace(-1,-10,10); % The equilibrium voltage of the AgLi particle in the surface defect, in unit mV; 0 corresponds to Li metal anode;

% nEn   =0;
for jLi = 1 : length(sgmLi)
    for jEn = 1 : length(sgmEn)
%         nEn = nEn+1;
%% Read data from files
% Relation of filenames and variable names: "sgm001" and "sgm-2" means sgm = 10^-2 = 0.01 mS/cm 
% Load potentials at different interfaces
    tmpStr  = ['sePot_Li',num2str(jLi),'En',num2str(jEn),'= csvread(''../rst/sgm_sePot_Li',num2str(jLi),'En',num2str(jEn),'.csv'',1,0);'];
    eval(tmpStr);
    tmpStr  = ['porePot_Li',num2str(jLi),'En',num2str(jEn),'= csvread(''../rst/sgm_porePot_Li',num2str(jLi),'En',num2str(jEn),'.csv'',1,0);'];
    eval(tmpStr);
    tmpStr  = ['agPot_Li',num2str(jLi),'En',num2str(jEn),'= csvread(''../rst/sgm_agPot_Li',num2str(jLi),'En',num2str(jEn),'.csv'',1,0);'];
    eval(tmpStr);
    tmpStr  = ['ccPot_Li',num2str(jLi),'En',num2str(jEn),'= csvread(''../rst/sgm_ccPot_Li',num2str(jLi),'En',num2str(jEn),'.csv'',1,0);'];
    eval(tmpStr);
%     tmpStr  = ['x1Pot_Li',num2str(jLi),'En',num2str(jEn),'= csvread(''../rst/sgm_z1Pot_Li',num2str(jLi),'En',num2str(jEn),'.csv'',1,0);'];
%     eval(tmpStr);
%     tmpStr  = ['x2Pot_Li',num2str(jLi),'En',num2str(jEn),'= csvread(''../rst/sgm_z2Pot_Li',num2str(jLi),'En',num2str(jEn),'.csv'',1,0);'];
%     eval(tmpStr);
%     tmpStr  = ['x3Pot_Li',num2str(jLi),'En',num2str(jEn),'= csvread(''../rst/sgm_z3Pot_Li',num2str(jLi),'En',num2str(jEn),'.csv'',1,0);'];
%     eval(tmpStr);

%% Data analysis
% Compute normal current density from the analytic Butler-Volmer relation:
 % i_n = i_exc * (exp(alpha*eta*F/RT)-exp(-alpha*eta*F/RT))
    tmpStr  = ['seEta_Li',num2str(jLi),'En',num2str(jEn),'  = - sePot_Li',num2str(jLi),'En',num2str(jEn),'(:,2) - sePot_Li',num2str(jLi),'En',num2str(jEn),'(:,3);'];
    eval(tmpStr);
    tmpStr  = ['seCrnt_Li',num2str(jLi),'En',num2str(jEn),' = iSE_exc * (exp(a0*F_RT*seEta_Li',num2str(jLi),'En',num2str(jEn),') - exp(-a0*F_RT*seEta_Li',num2str(jLi),'En',num2str(jEn),'));'];
    eval(tmpStr);
    tmpStr  = ['poreEta_Li',num2str(jLi),'En',num2str(jEn),'  = - porePot_Li',num2str(jLi),'En',num2str(jEn),'(:,2) - porePot_Li',num2str(jLi),'En',num2str(jEn),'(:,3);']; 
    eval(tmpStr);
    tmpStr  = ['poreCrnt_Li',num2str(jLi),'En',num2str(jEn),' = iSE_exc * (exp(a0*F_RT*poreEta_Li',num2str(jLi),'En',num2str(jEn),') - exp(-a0*F_RT*poreEta_Li',num2str(jLi),'En',num2str(jEn),'));'];
    eval(tmpStr);
    tmpStr  = ['agEta_Li',num2str(jLi),'En',num2str(jEn),'  = - - agPot_Li',num2str(jLi),'En',num2str(jEn),'(:,2) - agPot_Li',num2str(jLi),'En',num2str(jEn),'(:,3);']; 
    eval(tmpStr);
    tmpStr  = ['agCrnt_Li',num2str(jLi),'En',num2str(jEn),' = iAg_exc * (exp(a0*F_RT*agEta_Li',num2str(jLi),'En',num2str(jEn),') - exp(-a0*F_RT*agEta_Li',num2str(jLi),'En',num2str(jEn),'));'];
    eval(tmpStr);
    tmpStr  = ['ccEta_Li',num2str(jLi),'En',num2str(jEn),'  = - ccPot_Li',num2str(jLi),'En',num2str(jEn),'(:,2) - ccPot_Li',num2str(jLi),'En',num2str(jEn),'(:,3);'];
    eval(tmpStr);
    tmpStr  = ['ccCrnt_Li',num2str(jLi),'En',num2str(jEn),' = iCC_exc * (exp(a0*F_RT*ccEta_Li',num2str(jLi),'En',num2str(jEn),') - exp(-a0*F_RT*ccEta_Li',num2str(jLi),'En',num2str(jEn),'));'];
    eval(tmpStr);

% % Compute the total current ratio at each interfaces
%     tmpStr   =['seIave(jAg) = sum(seCrnt_Li',num2str(jLi),'En',num2str(jEn),')/length(seCrnt_Li',num2str(jLi),'En',num2str(jEn),');'];  eval(tmpStr);
%     tmpStr   =['ccIave(jAg) = sum(ccCrnt_Li',num2str(jLi),'En',num2str(jEn),')/length(ccCrnt_Li',num2str(jLi),'En',num2str(jEn),');'];  eval(tmpStr);
%     tmpStr   =['poreIave(jAg) = sum(poreCrnt_Li',num2str(jLi),'En',num2str(jEn),')/length(poreCrnt_Li',num2str(jLi),'En',num2str(jEn),');'];  eval(tmpStr);
%     tmpStr   =['agIave(jAg) = sum(agCrnt_Li',num2str(jLi),'En',num2str(jEn),')/length(agCrnt_Li',num2str(jLi),'En',num2str(jEn),');'];  eval(tmpStr);

    % Compute the total current ratio at each interfaces
    tmpStr   =['IrSE(jLi,jEn) = -mean(seCrnt_Li',num2str(jLi),'En',num2str(jEn),');'];  eval(tmpStr);
    tmpStr   =['IrCC(jLi,jEn) = -mean(ccCrnt_Li',num2str(jLi),'En',num2str(jEn),');'];  eval(tmpStr);
    tmpStr   =['IrPr(jLi,jEn) = -mean(poreCrnt_Li',num2str(jLi),'En',num2str(jEn),')*pi*dia/yBL;'];  eval(tmpStr);
    tmpStr   =['IrAg(jLi,jEn) = -mean(agCrnt_Li',num2str(jLi),'En',num2str(jEn),')*pi*dia/yBL;'];  eval(tmpStr);


    end

end

%% Plotting
% Plot potential at sgm= 0.8 mS/cm, ASR= 2 Ohm*cm^2 as a function of x coordinate at four different interfaces
ifg = ifg + 1;
figure(ifg)
sLi = repmat(sgmLi',1,length(sgmEn));
sEn = repmat(sgmEn,length(sgmLi),1);
contourf(sEn,sLi,IrCC);

sgmRt = sEn ./ sLi;
sgmRt_plt = reshape(sgmRt,16*16,1);
IrCC_plt = reshape(IrCC,16*16,1);
IrSE_plt = reshape(IrSE,16*16,1);

ifg = ifg + 1;
figure(ifg)
plot(sgmRt_plt,IrCC_plt,'o');

ifg = ifg + 1;
figure(ifg)
plot(sgmRt_plt,IrSE_plt,'o');

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
