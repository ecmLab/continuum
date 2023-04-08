%% This code is developed for the Li-Ag project, and the corresponding paper
% This version considered the Liniearized BV relation, 
% The intraparticle phase separation supressed
% Howard Tu, 06/27/21
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
La          = 5;                % The distance of the particles to the anode, in unit um
dia         = 0.1;            % The diameter of particles, in unit um
Np          = length(dia);           % the number of particles included
% 3. Material properties in the AgC interlayer
sgm_Li      = 0.5;                   % Ionic conductivity, in mS/cm
sgm_e       = 1;                  % Electronic conductivity, in unit mS/cm
c_e         = 1;                  % The relative electrons concentration, dimensionless, due to electron migration limitation
% 4. Material properties of the Ag particles
i_exc       = 1.3;                   % According to Chiku's paper (Microelectrode studies on charger transfer), unit mA/cm^2
R_Mt        = 1000/(i_exc*FF/RT);    % ASR of Li-metal/SE charge-transfer resistance, in unit Ohm*cm^2
R_Pt        = R_Mt/c_e;              % The electron-limited ARS for Void-SE charge-transfer resistance, in unit Ohm*cm^2
R_And       = R_Mt;                  % ASR for anode/interflayer charge-transfer resistance, in unit Ohm*cm^2
yp_end      = 0.9/(1-0.9);           % The terminal value of y in AgLi_y based on the solubility limit x=0.885
eng_den     = yp_end*FF/(Ag_mmol+yp_end*Li_mmol)*10/36;   % The energy density of AgLi10, in unit mAh/g
mas_den     = 1.0;                   % Mass density of AgLi10, in unit g/cm^3
rho         = yp_end*mas_den/(Ag_mmol+yp_end*Li_mmol); % Interstitial site density of AgLi9, in unit mol/cm^3 = Mass_den / Molar_Mass
% 5. Electrochemical parameters
i_ave       = 0.2;                 % Applied current density, in unit mA/cm^2
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
for ip = 1 : Np
    M11(ip,ip:ip+1) = La(ip)*[1,-1];
    M12(ip,ip:ip+1) = [1,-1];
    M43(ip,ip)      = La(ip);
end
M42 = [eye(Np),zeros(Np,1)];     % In unit 1
M31 =  sgm_Li*M42;          % In unit mS/cm
M33 = -sgm_e*M42;           % In unit mS/cm
M41 = M43 + sgm_Li*R_Pt*M12*10;  % In unit um
M51 = zeros(4,Np+1);  M51(1,1)    = -sgm_Li*R_And*10;  M51(2,Np+1) = 1;
M52 = zeros(4,Np+1);  M52(1,1)    = 1;
M53 = zeros(4,Np+1);  M53(3,Np+1) = 1;
M54 = zeros(4,Np+1);  M54(4,1)    = 1;
% Assemble the MM matrix, size of MM: (4*Np+4) x (4*Np+4)
MM  = [M11, M12, M00, M00; M00, M00, M11, M12; M31, M00, M33, M00; M41, M42, M43, M42; M51, M52, M53, M54]; 
% Defind Y in the matrix equation: MM*X = Y. Size of Y vector: (4*Np+4) x 1
% Entry 1:2*Np, 4*Np+3, 4*Np+4 are always zero; Entry 2*Np+1:3*Np and 4*Np+2 are constant;
% Entry 3*Np+1:4*Np are Li chemical potential in Ag particles; Entry 4*Np+1 is the Li chemical potential in the anode;
YY  = [zeros(1,2*Np),FF*i_ave/10^4*ones(1,Np),zeros(1,Np),0,FF*i_ave/(10^4*sgm_Li),0,0]';

% 2. The DD matrix in the equation: yc' = DD*m
DD = zeros(Np,Np+1);
for ip = 1 : Np
    DD(ip,ip:ip+1) = 6*sgm_Li/(rho*FF^2*dia(ip))*[1,-1]*10^5;       % In unit mol/J * um/s
end
          
% 3. Compute and load the Li chemical potential of AgLiy, y from 0 to yp_end, in unit J/mol
% LiAg_freEng_chemPot;
load('chemP_AgLi_db.mat');

%% Numerical solution with Forward Euler method
for it = 1:Nt-1
    chem_pot(:,it)   = chemPot_AgLi(chemP_AgLi_db,chemP_flag,yc(:,it)); % The Li chemical potential of all particles at current step, in unit J/mol                                            
%     YY(3*Np+1:4*Np)  = chem_pot(:,it);
%     YY(4*Np+1)       = mu_And(it);                                      % For Li metal anode, the Li chemical potential is always 0, in unit J/mol
    XX               = MM\YY;                                           % Solve the main equation of all 4*Np+4 unknowns
    m_cnt(:,it)      = XX(1:Np+1);
    n_cnt(:,it)      = XX(Np+2:2*(Np+1));
    b_cnt(:,it)      = XX(2*(Np+1)+1:3*(Np+1));
    c_cnt(:,it)      = XX(3*(Np+1)+1:4*(Np+1));
    ytmp             = dt*DD*m_cnt(:,it) + yc(:,it);                    % update the Li fraction at time it
    ytmp(ytmp>1)     = 1;
    yc(:,it+1)       = ytmp;
end

%% Plot
% 1. Plot potentials
nnp = 1000;
xxp = linspace(0,Lx,nnp)';
idx = zeros(nnp,1);
xrg = [0;La;Lx];
for ip = 1 : Np+1
    idx(xxp>=xrg(ip)&xxp<=xrg(ip+1)) = ip;
end
mu_ionic = zeros(nnp,1);
mu_elec  = zeros(nnp,1);
ifg = ifg + 1;
for it = 1:Nt-1
    for ip = 1 : nnp
        mu_ionic(ip) = m_cnt(idx(ip),it)*xxp(ip) + n_cnt(idx(ip),it);
        mu_elec(ip)  = b_cnt(idx(ip),it)*xxp(ip) + c_cnt(idx(ip),it);
    end
    figure(ifg)
    plot(xxp,mu_ionic/FF*1000,'-r', xxp, mu_elec/FF*1000,'-b', xxp, (mu_ionic+mu_elec)/FF*1000,'--g')
    pause(0.1)
end
% 1. Plot the equilibrium potential
yp=linspace(0,1,1000)';
eqmVlt_plt = chemPot_AgLi(eqmV_AgLi_db,chemP_flag,yp);
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
    plot(cell_DOD,chm_pot(ip,:)+Voc);
end
hold off
xlabel('Cell SOC (%)')
ylabel('Voltage (V)')

%% Define material functions
% The Li chemical potential of AgLi_y vs. Li content (y), in unit J/mol
% chemP_flag=1 means thermo-voltage, chemP_flag=2 means overshoot voltage
function chemPot = chemPot_AgLi(chemP_AgLi_db,chemP_flag,x)
   chemP_db = chemP_AgLi_db(:,chemP_flag+1);
   chemPot = interp1(chemP_AgLi_db(:,1),chemP_db,x);
end
