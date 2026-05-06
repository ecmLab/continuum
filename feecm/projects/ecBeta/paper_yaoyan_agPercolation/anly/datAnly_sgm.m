clc; clear; myFigureSetting;

ifg = 0;
%Geometric parameter
xBL = 10;      % width of the BL, in unit um
yBL = 20;      % thickness of the BL, in unit um
diaAg = 0.01;  % diameter of pore and Ag, in unit um
%Electrochemical parameters
F_RT  = 0.03868;  % The combined constant F/RT when T=300K, in unit 1/mV
a0    = 0.5;      % Reaction rate for the charge transfer reaction, set as symmetric anodic and cathodic reaction for now
% For LPS vs Li-metal system
iSE_exc = 1.3;        % The exchange current density, 13 is based on the ASR=2 Ohm, in unit mA/cm^2; LPS/Li metal interface from paper
iAg_exc = 10;  % The exchange current density, 13 is based on the ASR=2 Ohm, in unit mA/cm^2; LPS/Li metal interface from paper
iCC_exc = iSE_exc; % The exchange current density, 13 is based on the ASR=2 Ohm, in unit mA/cm^2; LPS/Li metal interface from paper
i0      = 0.68;     % applied current density, in unit mA/cm^2
% For change variables
sgmLi = 10.^(linspace(-1.1,1.1,12));  % The ionic conductivity, set as variable for parameter studies, in unit mS/cm
sgmEn = 172;  % The ionic conductivity, set as variable for parameter studies, in unit mS/cm
% sgmLi = 0.8;  % The ionic conductivity, set as variable for parameter studies, in unit mS/cm
% sgmEn = 10.^[0.01:0.2:3.01];  % The ionic conductivity, set as variable for parameter studies, in unit mS/cm
vAg   = -1; % The equilibrium voltage of the AgLi particle in the surface defect, in unit mV; 0 corresponds to Li metal anode;

% nEn   =0;
for jLi = 1 : length(sgmLi)
    %         nEn = nEn+1;
    %% Read data from files
    % Relation of filenames and variable names: "sgm001" and "sgm-2" means sgm = 10^-2 = 0.01 mS/cm
    % Load potentials at different interfaces
    tmpStr  = ['sePot_Li',num2str(jLi),'= csvread(''../rst/sgm_sePot_Li',num2str(jLi),'.csv'',1,0);'];
    eval(tmpStr);
    tmpStr  = ['agPot_Li',num2str(jLi),'= csvread(''../rst/sgm_agPot_Li',num2str(jLi),'.csv'',1,0);'];
    eval(tmpStr);
    tmpStr  = ['ccPot_Li',num2str(jLi),'= csvread(''../rst/sgm_ccPot_Li',num2str(jLi),'.csv'',1,0);'];
    eval(tmpStr);

    %% Data analysis
    % Compute normal current density from the analytic Butler-Volmer relation:
    % i_n = i_exc * (exp(alpha*eta*F/RT)-exp(-alpha*eta*F/RT))
    tmpStr  = ['seEta_Li',num2str(jLi),'  = - sePot_Li',num2str(jLi),'(:,2) - sePot_Li',num2str(jLi),'(:,3);'];
    eval(tmpStr);
    tmpStr  = ['seCrnt_Li',num2str(jLi),' = iSE_exc * (exp(a0*F_RT*seEta_Li',num2str(jLi),') - exp(-a0*F_RT*seEta_Li',num2str(jLi),'));'];
    eval(tmpStr);
    tmpStr  = ['agEta_Li',num2str(jLi),'  = - vAg - agPot_Li',num2str(jLi),'(:,2) - agPot_Li',num2str(jLi),'(:,3);'];
    eval(tmpStr);
    tmpStr  = ['agCrnt_Li',num2str(jLi),' = iAg_exc * (exp(a0*F_RT*agEta_Li',num2str(jLi),') - exp(-a0*F_RT*agEta_Li',num2str(jLi),'));'];
    eval(tmpStr);
    tmpStr  = ['ccEta_Li',num2str(jLi),'  = - ccPot_Li',num2str(jLi),'(:,2) - ccPot_Li',num2str(jLi),'(:,3);'];
    eval(tmpStr);
    tmpStr  = ['ccCrnt_Li',num2str(jLi),' = iCC_exc * (exp(a0*F_RT*ccEta_Li',num2str(jLi),') - exp(-a0*F_RT*ccEta_Li',num2str(jLi),'));'];
    eval(tmpStr);

    % Compute the total current ratio at each interfaces
    tmpStr   =['seIr(jLi) = -mean(seCrnt_Li',num2str(jLi),')/i0;'];  eval(tmpStr);
    tmpStr   =['ccIr(jLi) = -mean(ccCrnt_Li',num2str(jLi),')/i0;'];  eval(tmpStr);

end

%% Plotting
sgmRt = sgmEn ./ sgmLi;
ifg = ifg + 1;
figure(ifg)
semilogx(sgmRt,(1-seIr-ccIr)*100);

