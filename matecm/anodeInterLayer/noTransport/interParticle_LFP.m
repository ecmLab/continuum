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
% Parameters for the model
i_exc       = 8.5*10^(-4);           % Exchange current density, in unit mA/cm^2
i_ave       = 0.06*i_exc;            % Applied current density, in unit mA/cm^2
eng_den     = 170;                   % The energy density of LiFePO4, in unit mAh/g
mas_den     = 3.45;                  % Mass density of LiFePO4, in unit g/cm^3
rho         = 0.0228;                % Interstitial site density of LiFePO4
dia         = [0.01, 0.0175]';        % The diameter of particles, in unit um
Np          = length(dia);           % the number of particles included
TotT        = eng_den*mas_den*sum(dia.^3)/(6*i_ave*sum(dia.^2))*0.36; % The total discharging time depends on the applied current, in unit s
Nt          = 10000;                 % The total timesteps
dt          = TotT/Nt;               % the time step, in unit s
Voc         = 3.42;                  % Open circuit voltage plateau, in unit V
% Initiation
xc          = zeros(Np,Nt);          % Initialize Li site fraction of each particle
xc(:,1)     = 0.02*ones(Np,1);       % The initial Li fraction is 2%
eqm_vlt     = zeros(Np,Nt);          % The equilibrium voltage of all particles, in unit V
eqm_vlt(:,1)= eqmVlt_LPF(xc(:,1));  % The initial equilibrium voltage of all particles, in unit V
cell_vlt    = zeros(1,Nt);           % The cell voltage, in unit V
% For plot
ifg         = 0; 

%% Numerical solution with Forward Euler method
% Constants in the numerical equation
cnst1 = 6*i_exc*dt/(rho*RT)*10;         % in unit um/V
cnst2 = RT/FF*(i_ave/i_exc);            % in unit V
cnst3 = sum(dia.^2);                    % in unit um^2
% Forward Euler
for it = 2:Nt 
    cnst4            = sum(eqm_vlt(:,it-1) .* dia.^2);   % in unit V*um^2
    cell_vlt(it-1)   = Voc - cnst2 + cnst4/cnst3;        % The cell voltage at time it-1, in unit V
    xc(:,it)         = cnst1*(cnst2 - cnst4/cnst3 + eqm_vlt(:,it-1))./dia + xc(:,it-1);  % update the Li fraction at time it
    eqm_vlt(:,it)    = eqmVlt_LPF(xc(:,it));             % The equilibrium voltage of all particles at time it, in unit V
end
cell_vlt(:,Nt)       = Voc - cnst2 + sum(eqm_vlt(:,Nt) .* dia.^2)/cnst3; % The Last cell voltage, in unit V

%% Plot
% 1. Plot the free energy
xp=linspace(0,1,1000);
freEng_plt = freEng_LPF(xp,FF);
ifg = ifg + 1;
figure(ifg)
plot(xp,freEng_plt);
xlabel('x in Li_xFePO_4')
ylabel('Gibbs free energy w.r.t LiFePO_4 (J/mol)')
% 2. Plot the equilibrium potential
eqmVlt_plt = eqmVlt_LPF(xp);
ifg = ifg + 1;
figure(ifg)
plot(xp,eqmVlt_plt);
xlabel('x in Li_xFePO_4')
ylabel('Equilibrium voltage (V)')
% 3. Plot DOD of each particle
cell_DOD = linspace(0,100,Nt)';
ifg = ifg + 1;
figure(ifg)
for ip = 1 : Np
    plot(cell_DOD,xc(ip,:)*100);
    hold on
end
hold off
xlabel('Cell DOD (%)')
ylabel('Particle DOD (%)')
% 4. Plot Voltage of the cell and each particle
ifg = ifg + 1;
figure(ifg)
plot(cell_DOD, cell_vlt);
for ip = 1 : Np
    hold on
    plot(cell_DOD,eqm_vlt(ip,:)+Voc);
end
hold off
xlabel('Cell DOD (%)')
ylabel('Voltage (V)')

%% Define material functions
% 1. The equilibrium voltage of Li_xFePO4 vs. Li content (x), in unit V
% The eqmVolt = -chemical_potential/FF
function eqmVlt = eqmVlt_LPF(x)
    eqmVlt = (5*(1.02-2.04*x).^51 - 2.925275*x.^2 + 6.37507*x - 2.558325)*0.01;
end

% 2. The Gibbs free energy of Li_xFePO4 vs. Li content (x), in unit J/mol
% The freEng = integrate(chemical_potential)dx
function freEng = freEng_LPF(x,FF)
    freEng = -0.01*FF*(-5/(52*2.04)*(1.02-2.04*x).^52 - 2.925275/3*x.^3 + 6.37507/2*x.^2 - 2.558325*x)-461.0773;
end