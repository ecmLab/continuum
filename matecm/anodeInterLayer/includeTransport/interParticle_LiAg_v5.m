%% This code is developed for the Li-Ag project, for SE-MIEC
% This version considered the Liniearized BV relation, 
% The intraparticle phase separation supressed
% Howard Tu, 07/16/21
clc;clear;
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontSize',16);
ifg         = 0;

%% Parameters
% 1. Constants
FF          = 96485.33;              % Faraday constant, in unit s*A/mol
RT          = 8.314*300;             % RT constant, in unit J/mol
Ag_mmol     = 107.868;               % Molar mass of Silver, in unit g/mol
Li_mmol     = 6.948;                 % Molar mass of Lithium, in unit g/mol
% 2. Microstructural parameters
Lx          = 10;                    % Length of the interlayer, in unit um
La          = [1,2,3,4,5]';                     % The distance of the particles to the anode, in unit um
% La          = 0.5;                     % The distance of the particles to the anode, in unit um
dia         = [0.1,0.2,0.3,0.4,0.5]';           % The diameter of particles, in unit um
% dia         = 0.1;           % The diameter of particles, in unit um
Np          = length(dia);           % the number of particles included
% 3. Material properties in the AgC interlayer
sgm_Li      = 0.5;                   % Ionic conductivity, in mS/cm
c_e         = 1;                  % The relative electrons concentration, dimensionless, due to electron migration limitation
% 4. Material properties of the Ag particles
i_exc       = 1.3;                   % According to Chiku's paper (Microelectrode studies on charger transfer), unit mA/cm^2
R_Mt        = 1000/(i_exc*FF/RT);    % ASR of Li-metal/SE charge-transfer resistance, in unit Ohm*cm^2
R_Pt        = R_Mt/c_e;              % The electron-limited ARS for Void-SE charge-transfer resistance, in unit Ohm*cm^2
R_And       = R_Mt;                  % ASR for anode/interflayer charge-transfer resistance, in unit Ohm*cm^2
R_Ctd       = R_Mt;                  % ASR for cathode/interflayer charge-transfer resistance, in unit Ohm*cm^2
yp_end      = 0.9/(1-0.9);           % The terminal value of y in AgLi_y based on the solubility limit x=0.885
eng_den     = yp_end*FF/(Ag_mmol+yp_end*Li_mmol)*10/36;   % The energy density of AgLi10, in unit mAh/g
mas_den     = 1.0;                   % Mass density of AgLi10, in unit g/cm^3
rho         = yp_end*mas_den/(Ag_mmol+yp_end*Li_mmol); % Interstitial site density of AgLi9, in unit mol/cm^3 = Mass_den / Molar_Mass
% 5. Electrochemical parameters
i_ave       = 0.02;                 % Applied current density, in unit mA/cm^2
chemP_flag  = 1;                    % The chemical potential used, flag=1 means thermo-voltage, flag=2 means overshoot voltage
% 6. Temperal parameters
Nt          = 1000;                 % The total timesteps
TotT        = 10;                    % The total discharging time depends on the applied current, in unit s
dt          = TotT/Nt;               % the time step, in unit s
% 7. Initialization
yc          = zeros(Np,Nt);          % Initialize Li site fraction of each particle
yc(:,1)     = 0.5*ones(Np,1);        % The initial Li fraction is 2%
chem_pot    = zeros(Np,Nt);          % The Li chemical potential of all particles, in unit J/mol
mu_And      = zeros(1,Nt);           % The Li chemical potential in anode, in unit J/mol
mu_Ctd      = zeros(1,Nt);           % The Li chemical potential in cathode, in unit J/mol
m_cnt       = zeros(Np+1,Nt);        % The first constant in the ionic potential equation: u_M+ = m_cnt*x + n_cnt. In unit J/mol/um
n_cnt       = zeros(Np+1,Nt);        % The second constant in the ionic potential equation: u_M+ = m_cnt*x + n_cnt. In unit J/mol
V0          = zeros(1,Nt);           % The voltage of the cell. In unit J/mol

%% Compute all the matrix constant for solving equations
% 1. Components of MM in the matrix equation: MM*X = Y
M11 = zeros(Np,Np+1);  % In unit um
M12 = zeros(Np,Np+1);  % In unit 1
M13 = zeros(Np,1);
M21 = zeros(Np,Np+1);  % In unit um
for ip = 1 : Np
    M11(ip,ip:ip+1) = La(ip)*[1,-1];
    M12(ip,ip:ip+1) = [1,-1];
    M21(ip,ip)      = La(ip);
end
M21 = M21 + 10*sgm_Li*R_Pt*M12;   % In unit um
M22 = [eye(Np),zeros(Np,1)];      % In unit 1
M23 = ones(Np,1);

M31 = zeros(3,Np+1);  M31(1,1)    = sgm_Li*R_And*10;   M31(2,Np+1) = Lx+sgm_Li*R_Ctd*10;   M31(3,Np+1) = Lx;
M32 = zeros(3,Np+1);  M32(1,1)    = -1;                M32(2,Np+1) = 1;                    M32(3,Np+1) = 1;
M33 = zeros(3,1);     M33(1)      = -1;                M33(2)      = 1;
% Assemble the MM matrix, size of MM: (2*Np+3) x (2*Np+3)
MM  = [M11, M12, M13; M21, M22, M23; M31, M32, M33]; 
% Defind Y in the matrix equation: MM*X = Y. Size of Y vector: (2*Np+3) x 1
% Entry 1:Np and 2*Np+3 are always zero;
% Entry Np+1:2*Np are Li chemical potential in Ag particles; 
% Entry 2*Np+1 is the Li chemical potential in the anode; Entry 2*Np+2 is the Li chemical potential in the cathode;
YY  = [zeros(1,Np),chem_pot(:,1)',-mu_And(1)-FF*i_ave*R_And/10^3,mu_Ctd(1),0]';

% 2. The DD matrix in the equation: yc' = DD*m
DD  = zeros(Np,Np+1);
cnt = -10^5*sgm_Li/(rho*FF^2);
for ip = 1 : Np
    DD(ip,ip:ip+1) = cnt*[1,-1]/dia(ip);       % In unit mol/J * um/s
end

% 3. Compute and load the Li chemical potential of AgLiy, y from 0 to yp_end, in unit J/mol
% LiAg_freEng_chemPot;
load('chemP_AgLi_db.mat');

%% Numerical solution with Forward Euler method
for it = 1:Nt-1
    chem_pot(:,it)   = chemPot_AgLi(chemP_AgLi_db,chemP_flag,yc(:,it)); % The Li chemical potential of all particles at current step, in unit J/mol                                            
    YY(Np+1:2*Np)    = chem_pot(:,it);
    YY(2*Np+1)       = -mu_And(it)-FF*i_ave*R_And/10^3;                 % For Li metal anode, the Li chemical potential is always 0, in unit J/mol
    YY(2*Np+2)       = mu_Ctd(it);                                      % For Li metal cathode, the Li chemical potential is always 0, in unit J/mol 
    XX               = MM\YY;                                           % Solve the main equation of all 4*Np+4 unknowns
    m_cnt(:,it)      = XX(1:Np+1);
    n_cnt(:,it)      = XX(Np+2:2*Np+2);
    V0(it)           = XX(end);
    ytmp             = dt*DD*m_cnt(:,it) + yc(:,it);                    % update the Li fraction at time it
    ytmp(ytmp>1)     = 1;
    yc(:,it+1)       = ytmp;
end

%% Data analysis
% Deposition current density at the interfaces and in the particles. i=1~Anode, i=2:Np+2~Particles, i=Np+3~Cathode
depI            = zeros(Np+2,Nt-1);      
depI(1,:)       = i_ave + sgm_Li/FF*m_cnt(1,1:Nt-1)*10^4;   % Current density comsumed at the anode interface
depI(2:end-1,:) = -sgm_Li/FF*(m_cnt(1:end-1,1:Nt-1)-m_cnt(2:end,1:Nt-1))*10^4;    % Current density comsumed at the cathode interface
depI(end,:)     = -sgm_Li/FF*m_cnt(end,1:Nt-1)*10^4;    % Current density comsumed at the cathode interface

%% Plot
% 1. Plot potentials and current at the given timestep
tmStp = [1,10];
nnp = 1000;
xxp = linspace(0,Lx,nnp)';
idx = zeros(nnp,1);
xrg = [0;La;Lx];
for ip = 1 : Np+1
    idx(xxp>=xrg(ip)&xxp<=xrg(ip+1)) = ip;
end
mu_ionic = zeros(nnp,1);
it_ionic = zeros(nnp,1);
ifg = ifg + 1;
for it = 1:length(tmStp)
    for ip = 1 : nnp
        mu_ionic(ip) = m_cnt(idx(ip),tmStp(it))*xxp(ip) + n_cnt(idx(ip),tmStp(it));
        it_ionic(ip) = -m_cnt(idx(ip),tmStp(it))*sgm_Li/FF*10^4;
    end
    figure(ifg)
    plot(xxp-Lx,mu_ionic/FF*10^3,'-r', xxp-Lx, V0(tmStp(it))/FF*10^3,'-b', xxp-Lx, (mu_ionic+V0(tmStp(it)))/FF*10^3,'--g')
    xlabel('Distance to the CC (\mum)');
    ylabel('Voltage (\muV)');
    legend('Ionic','Electronic','Plating');
    figure(ifg+1)
    plot(xxp-Lx,it_ionic/i_ave,'-r')
    xlabel('Distance to the CC (\mum)');
    ylabel('Current density (mA/cm^2)');
    
%     pause;
end

% 2. Plot the current percentage at the interface and at the particles as a function of time
ttp = linspace(dt,TotT,Nt-1)';
ifg = ifg + 2;
figure(ifg)
for ii = 1 : length(depI(:,1))
    plot(ttp,depI(ii,:)/i_ave);
    hold on
end
hold off
xlabel('Charging time (s)');
ylabel('Current density rate');
legend('Anode','Particle1','Particle2','Cathode');
    
% 3. Plot the voltage as a function of time
ifg = ifg + 1;
figure(ifg)
plot(ttp,-V0(1:Nt-1)/FF*10^3);
xlabel('Charging time (s)');
ylabel('Voltage (mV)');

%% Define material functions
% The Li chemical potential of AgLi_y vs. Li content (y), in unit J/mol
% chemP_flag=1 means thermo-voltage, chemP_flag=2 means overshoot voltage
function chemPot = chemPot_AgLi(chemP_AgLi_db,chemP_flag,x)
   chemP_db = chemP_AgLi_db(:,chemP_flag+1);
   chemPot = interp1(chemP_AgLi_db(:,1),chemP_db,x);
end
