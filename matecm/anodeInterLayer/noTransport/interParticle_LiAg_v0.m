%% This code is developed for the Li-Ag project, and the corresponding paper
% This version considered the Liniearized BV relation, 
% The intraparticle phase separation supressed
% Howard Tu, 04/27/21
clc;clear;
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontSize',16);

%% Parameters
% Constants
FF          = 96485.33;              % Faraday constant, in unit s*A/mol
RT          = 8.314*300;             % RT constant, in unit J/mol
Ag_mmol     = 107.868;               % Molar mass of Silver, in unit g/mol
Li_mmol     = 6.948;                 % Molar mass of Lithium, in unit g/mol
% Parameters for the model
i_exc       = 8.5*10^(-4);           % Exchange current density, in unit mA/cm^2
i_ave       = 0.06*i_exc;            % Applied current density, in unit mA/cm^2
yp_end      = 0.9/(1-0.9);       % The terminal value of y in AgLi_y based on the solubility limit x=0.885
eng_den     = yp_end*FF/(Ag_mmol+yp_end*Li_mmol)*10/36;   % The energy density of AgLi10, in unit mAh/g
mas_den     = 1.0;                  % Mass density of AgLi10, in unit g/cm^3
rho         = yp_end*mas_den/(Ag_mmol+yp_end*Li_mmol); % Interstitial site density of AgLi9, in unit mol/cm^3 = Mass_den / Molar_Mass
dia         = [0.01,0.02]';        % The diameter of particles, in unit um
% dia         = [0.01,0.02, 0.04,0.06,0.08]';        % The diameter of particles, in unit um
Np          = length(dia);           % the number of particles included
TotT        = eng_den*mas_den*sum(dia.^3)/(6*i_ave*sum(dia.^2))*0.36; % The total discharging time depends on the applied current, in unit s
Nt          = 10000;                 % The total timesteps
dt          = TotT/Nt;               % the time step, in unit s
Voc         = 0;                  % Open circuit voltage plateau, in unit V
ifg         = 0;
% The equilibrium voltage used in the calculation, flag=1 means thermo-voltage, flag=2 means overshoot voltage
eqmV_flag   = 1;                  
% Compute and load the equilibrium voltage of AgLiy, y from 0 to yp_end, in unit V
LiAg_freEng_eqmVlt;
load('eqmV_AgLi_db.mat');
% Initiation
yc          = zeros(Np,Nt);          % Initialize Li site fraction of each particle
yc(:,1)     = 0.02*ones(Np,1);       % The initial Li fraction is 2%
eqm_vlt     = zeros(Np,Nt);          % The equilibrium voltage of all particles, in unit V
eqm_vlt(:,1)= eqmVlt_AgLi(eqmV_AgLi_db,eqmV_flag,yc(:,1));  % The initial equilibrium voltage of all particles, in unit V
cell_vlt    = zeros(1,Nt);           % The cell voltage, in unit V

%% Numerical solution with Forward Euler method
% Constants in the numerical equation
cnst1 = 6*i_exc*dt/(rho*RT)*10;         % in unit um/V
cnst2 = RT/FF*(i_ave/i_exc);            % in unit V
cnst3 = sum(dia.^2);                    % in unit um^2
% Forward Euler
for it = 2:Nt 
    cnst4            = sum(eqm_vlt(:,it-1) .* dia.^2);   % in unit V*um^2
    cell_vlt(it-1)   = Voc - cnst2 + cnst4/cnst3;        % The cell voltage at time it-1, in unit V
    ytmp             = cnst1*(cnst2 - cnst4/cnst3 + eqm_vlt(:,it-1))./dia + yc(:,it-1);  % update the Li fraction at time it
    ytmp(ytmp>1)     = 1;
    yc(:,it)         = ytmp;
    eqm_vlt(:,it)    = eqmVlt_AgLi(eqmV_AgLi_db,eqmV_flag,yc(:,it));  % The equilibrium voltage of all particles at time it, in unit V
end
cell_vlt(:,Nt)       = Voc - cnst2 + sum(eqm_vlt(:,Nt) .* dia.^2)/cnst3; % The Last cell voltage, in unit V

%% Plot
% 1. Plot the equilibrium potential
yp=linspace(0,1,1000)';
eqmVlt_plt = eqmVlt_AgLi(eqmV_AgLi_db,eqmV_flag,yp);
ifg = ifg + 1;
figure(ifg)
plot(yp*100,eqmVlt_plt);
xlabel('Single Particle SOC (%)')
ylabel('Equilibrium voltage (V)')
% 2. Plot DOD of each particle
cell_DOD = linspace(0,100,Nt)';
ifg = ifg + 1;
figure(ifg)
for ip = 1 : Np
    plot(cell_DOD,yc(ip,:)*100);
    hold on
end
hold off
xlabel('Cell SOC (%)')
ylabel('Particle SOC (%)')
axis([0,100,0,105]);
% 3. Plot Voltage of the cell and each particle
ifg = ifg + 1;
figure(ifg)
plot(cell_DOD, cell_vlt);
for ip = 1 : Np
    hold on
    plot(cell_DOD,eqm_vlt(ip,:)+Voc);
end
hold off
xlabel('Cell SOC (%)')
ylabel('Voltage (V)')

%% Define material functions
% The equilibrium voltage of AgLi_y vs. Li content (y), in unit V
% eqmV_flag=1 means thermo-voltage, eqmV_flag=2 means overshoot voltage
function eqmVlt = eqmVlt_AgLi(eqmV_AgLi_db,eqmV_flag,x)
   eqmV_db = eqmV_AgLi_db(:,eqmV_flag+1);
   eqmVlt = interp1(eqmV_AgLi_db(:,1),eqmV_db,x);
end
