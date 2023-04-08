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
Hx          = 0.2;                   % The width of the interlayer, in unit um
La          = 1;                     % The distance of the particles to the anode, in unit um
dia         = 0.1;                   % The diameter of particles, in unit um
Np          = length(dia);           % the number of particles included
% 3. Material properties in the AgC interlayer
sgm_Li      = 0.5;                   % Ionic conductivity, in mS/cm
sgm_e       = 5;                  % Electronic conductivity, in unit mS/cm
c_e         = 0.1;                  % The relative electrons concentration, dimensionless, due to electron migration limitation
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
chemP_flag  = 1; % The chemical potential used, flag=1 means thermo-voltage, flag=2 means overshoot voltage
% 6. Temperal parameters
Nt          = 1000;                 % The total timesteps
TotT        = 10;                    % The total discharging time depends on the applied current, in unit s
dt          = TotT/Nt;               % the time step, in unit s
% 7. Initialization
yc          = zeros(Np,Nt);          % Initialize Li site fraction of each particle
yc(:,1)     = 0.02*ones(Np,1);       % The initial Li fraction is 2%
chem_pot    = zeros(Np,Nt);          % The Li chemical potential of all particles, in unit J/mol
mu_And      = zeros(1,Nt);           % The Li chemical potential in anode, in unit J/mol
mu_Ctd      = zeros(1,Nt);           % The Li chemical potential in cathode, in unit J/mol
m_cnt       = zeros(Np+1,Nt);        % The first constant in the ionic potential equation: u_M+ = m_cnt*x + n_cnt. In unit J/mol/um
n_cnt       = zeros(Np+1,Nt);        % The second constant in the ionic potential equation: u_M+ = m_cnt*x + n_cnt. In unit J/mol
b_cnt       = zeros(Np+1,Nt);        % The first constant in the electronic potential equation: u_e- = b_cnt*x + c_cnt. In unit J/mol/um
c_cnt       = zeros(Np+1,Nt);        % The second constant in the electronic potential equation: u_e- = b_cnt*x + c_cnt. In unit J/mol

%% Compute all the matrix constant for solving equations
% 1. Components of MM in the matrix equation: MM*X = Y
M00 = zeros(Np,Np+1);
M11 = zeros(Np,Np+1);  % In unit um
M12 = zeros(Np,Np+1);  % In unit 1
M43 = zeros(Np,Np+1);  % In unit um
Mtm = zeros(Np,Np+1);  % Temperary matrix for calculation
for ip = 1 : Np
    M11(ip,ip:ip+1) = La(ip)*[1,-1];
    M12(ip,ip:ip+1) = [1,-1];
    M43(ip,ip)      = La(ip);
    Mtm(ip,ip:ip+1) = [1,-1]/(1+2*dia(ip)/Hx);
end
M42 = [eye(Np),zeros(Np,1)];     % In unit 1
M41 = M43 + 5*sgm_Li*R_Pt*Mtm;   % In unit um
M31 = sgm_Li*eye(Np+1);          % In unit mS/cm
M32 = zeros(Np+1,Np+1);
M33 = -sgm_e*eye(Np+1);           % In unit mS/cm
M51 = zeros(3,Np+1);  M51(1,1)    = sgm_Li*R_And*10;   M51(2,Np+1) = Lx+sgm_Li*R_Ctd*10;
M52 = zeros(3,Np+1);  M52(1,1)    = -1;                M52(2,Np+1) = 1;
M53 = zeros(3,Np+1);  M53(3,Np+1) = Lx;
M54 = zeros(3,Np+1);  M54(1,1)    = -1;                M54(3,Np+1) = 1;
% Assemble the MM matrix, size of MM: (4*Np+4) x (4*Np+4)
MM  = [M11, M12, M00, M00; M00, M00, M11, M12; M31, M32, M33, M32; M41, M42, M43, M42; M51, M52, M53, M54]; 
% Defind Y in the matrix equation: MM*X = Y. Size of Y vector: (4*Np+4) x 1
% Entry 1:2*Np and 4*Np+4 are always zero; Entry 2*Np+1:3*Np+1 are constant;
% Entry 3*Np+2:4*Np+1 are Li chemical potential in Ag particles; 
% Entry 4*Np+2 is the Li chemical potential in the anode; Entry 4*Np+3 is the Li chemical potential in the cathode;
YY  = [zeros(1,2*Np),-FF*i_ave/10^4*ones(1,Np+1),chem_pot(:,1)',-mu_And(1)-FF*i_ave*R_And/10^3,mu_Ctd(1),0]';

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
%     YY(3*Np+2:4*Np+1)= chem_pot(:,it);
    YY(4*Np+2)       = -mu_And(it)-FF*i_ave*R_And/10^3;                 % For Li metal anode, the Li chemical potential is always 0, in unit J/mol
    YY(4*Np+3)       = mu_Ctd(it);                                      % For Li metal cathode, the Li chemical potential is always 0, in unit J/mol 
    XX               = MM\YY;                                           % Solve the main equation of all 4*Np+4 unknowns
    m_cnt(:,it)      = XX(1:Np+1);
    n_cnt(:,it)      = XX(Np+2:2*(Np+1));
    b_cnt(:,it)      = XX(2*(Np+1)+1:3*(Np+1));
    c_cnt(:,it)      = XX(3*(Np+1)+1:4*(Np+1));
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
tmStp = [1,2];
nnp = 1000;
xxp = linspace(0,Lx,nnp)';
idx = zeros(nnp,1);
xrg = [0;La;Lx];
for ip = 1 : Np+1
    idx(xxp>=xrg(ip)&xxp<=xrg(ip+1)) = ip;
end
mu_ionic = zeros(nnp,1);
mu_elec  = zeros(nnp,1);
it_ionic = zeros(nnp,1);
it_elec  = zeros(nnp,1);
ifg = ifg + 1;
for it = 1:length(tmStp)
    for ip = 1 : nnp
        mu_ionic(ip) = m_cnt(idx(ip),tmStp(it))*xxp(ip) + n_cnt(idx(ip),tmStp(it));
        mu_elec(ip)  = b_cnt(idx(ip),tmStp(it))*xxp(ip) + c_cnt(idx(ip),tmStp(it));
        it_ionic(ip) = -m_cnt(idx(ip),tmStp(it))*sgm_Li/FF*10^4;
        it_elec(ip)  = b_cnt(idx(ip),tmStp(it))*sgm_e/FF*10^4;
    end
    figure(ifg)
    plot(xxp-Lx,mu_ionic/FF*1000,'-r', xxp-Lx, mu_elec/FF*1000,'-b', xxp-Lx, (mu_ionic+mu_elec)/FF*1000,'--g')
    xlabel('Distance to the CC (\mum)');
    ylabel('Voltage (mV)');
    legend('Ionic','Electronic','Plating');
    figure(ifg+1)
    plot(xxp-Lx,it_ionic/i_ave,'-r', xxp-Lx, it_elec/i_ave,'-b')
    xlabel('Distance to the CC (\mum)');
    ylabel('Current density (mA/cm^2)');
    legend('Ionic','Electronic');
    
    pause;
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
legend('Anode','Particle','Cathode');
    
%% Define material functions
% The Li chemical potential of AgLi_y vs. Li content (y), in unit J/mol
% chemP_flag=1 means thermo-voltage, chemP_flag=2 means overshoot voltage
function chemPot = chemPot_AgLi(chemP_AgLi_db,chemP_flag,x)
   chemP_db = chemP_AgLi_db(:,chemP_flag+1);
   chemPot = interp1(chemP_AgLi_db(:,1),chemP_db,x);
end
