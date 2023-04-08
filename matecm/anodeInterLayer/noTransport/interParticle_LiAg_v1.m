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
Li_vmol     = 12.97;                 % Molar volume of Lithium, in unit cm^3/mol
% 2. Microstructural parameters
Lx          = 10;                    % Length of the interlayer, in unit um
Hx          = 0.2;                   % cross-section width, in unit um
% La          = [1,2,3,4,5]';                     % The distance of the particles to the anode, in unit um
La          = 0.5;                     % The distance of the particles to the anode, in unit um
% dia         = [0.1,0.2,0.3,0.4,0.5]';           % The diameter of particles, in unit um
dia         = [0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15];      % The diameter of particles, in unit um
Np          = length(dia);                 % the number of particles included
rVol        = pi/6*sum(dia.^3)/(Lx*Hx^2);  % The approximated volume ratio of Ag paritcles
% 3. Material properties in the AgC interlayer
sgm_Li      = 0.5;                   % Ionic conductivity, in mS/cm
c_p         = 1;
c_SE        = 0.1;
% 4. Material properties of the Ag particles
i_exc       = 0.13;                   % According to Chiku's paper (Microelectrode studies on charger transfer), unit mA/cm^2
R_Mt        = 1000/(i_exc*FF/RT);    % ASR of Li-metal/SE charge-transfer resistance, in unit Ohm*cm^2
R_Pt        = R_Mt/c_p;              % The electron-limited ARS for Void-SE charge-transfer resistance, in unit Ohm*cm^2
R_SE        = R_Mt/c_SE;                  % ASR for the SE/interflayer charge-transfer resistance, in unit Ohm*cm^2
R_CC        = R_Mt;             % ASR for the Current Collector/interflayer charge-transfer resistance, in unit Ohm*cm^2
yp_end      = 0.99/(1-0.99);         % The terminal value of y in AgLi_y based on the solubility limit x=0.885
eng_den     = yp_end*FF/(Ag_mmol+yp_end*Li_mmol)*10/36;   % The energy density of AgLi10, in unit mAh/g
mas_den     = massDensity_AgLi(yp_end);                   % Mass density of AgLi10, in unit g/cm^3
rho         = yp_end*mas_den/(Ag_mmol+yp_end*Li_mmol); % Interstitial site density of AgLi9, in unit mol/cm^3 = Mass_den / Molar_Mass
% 5. Electrochemical parameters
i_ave       = 2;                 % Applied current density, in unit mA/cm^2
chemP_flag  = 1;                    % The chemical potential used, flag=1 means thermo-voltage, flag=2 means overshoot voltage
% 6. Temperal parameters
Nt          = 10000;                 % The total timesteps
TotT        = 500;                    % The total discharging time depends on the applied current, in unit s
dt          = TotT/Nt;               % the time step, in unit s
% 7. Initialization
yc          = zeros(Np,Nt);          % Initialize Li site fraction of each particle, and at the interfaces
yc(:,1)     = 0.01*ones(Np,1);        % The initial Li fraction in Ag particles is 2%
chem_pot    = zeros(Np,Nt);          % The Li chemical potential of all particles, in unit J/mol
mu_SE       = zeros(1,Nt);           % The Li chemical potential at the SE/AgC interface, in unit J/mol
mu_CC       = zeros(1,Nt);           % The Li chemical potential at the CC/AgC interface, in unit J/mol
depI        = zeros(Np,Nt);          % The deposition current in the particles, in unit mA/cm^2
depSE       = zeros(1,Nt);           % The deposition current at SE/AgC interface, in unit mA/cm^2
depCC       = zeros(1,Nt);           % The deposition current at CC/AgC interface, in unit mA/cm^2
LiSE        = zeros(1,Nt);           % The thickness of deposited Li at the SE/AgC interface, in unit um
LiCC        = zeros(1,Nt);           % The thickness of deposited Li at the CC/AgC interface, in unit um
V0          = zeros(1,Nt);           % The voltage of the cell. In unit J/mol

%% Compute all the matrix constant for solving equations
% 1. Components of MM in the matrix equation: MM*X = Y
M11            = -FF*R_Pt/10^3*eye(Np+2);  % In unit J/mol/(mA/cm^2)
M11(Np+1,Np+1) = -FF*R_SE/10^3;
M11(Np+2,Np+2) = -FF*R_CC/10^3;
M12            = ones(Np+2,1);             % in unit 1
M21            = ones(1,Np+2);             % In unit 1
for ip = 1 : Np
    M21(ip)    = pi*(dia(ip)/Hx)^2;
end
M22 = 0;
% Assemble the MM matrix, size of MM: (Np+3) x (Np+3)
MM  = [M11, M12; M21, M22];
% Defind Y in the matrix equation: MM*X = Y. Size of Y vector: (Np+3) x 1
% Entry 1:Np are Li chemical potential in Ag particles; 
% Entry Np+1 is the Li chemical potential at SE/AgC; Entry Np+2 is the Li chemical potential atthe CC/AgC;
% Entry Np+3 is the cell voltage
YY  = [chem_pot(:,1)',mu_SE(1),mu_CC(1),i_ave]';
XX  = zeros(Np+3,1);

% 2. The DD matrix in the equation: yc' = DD*m
DD  = zeros(Np,Np);
cnt = 10*6/(rho*FF);
for ip = 1 : Np
    DD(ip,ip) = cnt/dia(ip);       % In unit 1/(mA/cm^2) * (1/s)
end

% 3. Compute and load the Li chemical potential of AgLiy, y from 0 to yp_end, in unit J/mol
% LiAg_freEng_chemPot;
load('chemP_AgLi_db.mat');

%% Numerical solution with Forward Euler method
for it = 1:Nt-1
%     if it == Nt-10
%         it
%     end
% Step 1: Try the first calculation
    for ip = 1:Np
        if yc(ip,it) < 1
            chem_pot(ip,it)   = chemPot_AgLi(chemP_AgLi_db,chemP_flag,yc(ip,it)); % The Li chemical potential in Ag particles at current step, in unit J/mol
        else
            chem_pot(ip,it)   = chemP_AgLi_db(end,chemP_flag+1);         % If Li content is 100%, Li chemical potential equals to Li metal chemical potential
        end
    end
    YY(1:Np)         = chem_pot(:,it);
    YY(Np+1)         = mu_SE(it);                                       % For Li metal anode, the Li chemical potential is always 0, in unit J/mol
    YY(Np+2)         = mu_CC(it);                                       % For Li metal cathode, the Li chemical potential is always 0, in unit J/mol 
    XX               = MM\YY;                                           % Solve the main equation of all Np+3 unknowns

% Step 2: Compute the accumulate Li content in the Ag particles and at the interfaces
    ytmp           = dt*DD*XX(1:Np) + yc(:,it);                      % update the Li fraction at time it
    SEtmp          = 10*dt*Li_vmol/FF*XX(Np+1) + LiSE(it);           % update the Li thickness at the SE/AgC at time it
    CCtmp          = 10*dt*Li_vmol/FF*XX(Np+2) + LiCC(it);           % update the Li thickness at the CC/AgC at time it

% Step 3: Impose constraints: 
% (1) Li content in particles cannot be less than zero; 
% (2) Li thickness at interfaces cannot be less than zero
    [tmp,yc0]      = find(ytmp<0);
%     [tmp,yc1]      = find(ytmp>1);
    [tmp,yct]      = find([SEtmp,CCtmp]<0);
    if (length(yc0)+length(yct)) > 0
       MMt         = MM;
       YYt         = YY;
       for iy = 1:length(yc0)
          MMt(yc0(iy),yc0(iy)) = 1;
          MMt(yc0(iy),Np+3)    = 0;
          YYt(yc0(iy))         = 0;
       end
%        for iy = 1:length(yc1)
%           YYt(yc1(iy))         = 0;
%        end
       for iy = 1:length(yct)
          MMt(yct(iy)+Np,yct(iy)+Np) = 1;
          MMt(yct(iy)+Np,Np+3)    = 0;
          YYt(yct(iy)+Np)         = 0;
       end
       XX               = MMt\YYt;                  % Solve the main equation again under constrain
    end
    
% Step 4: Compute the deposition rate and the cell voltage at it+1
    depI(:,it+1)     = XX(1:Np);
    depSE(it+1)      = XX(Np+1);
    depCC(it+1)      = XX(Np+2);
    V0(it+1)         = XX(Np+3); 
    yc(:,it+1)       = dt*DD*depI(:,it+1) + yc(:,it);
    LiSE(it+1)       = 10*dt*Li_vmol/FF*XX(Np+1) + LiSE(it);           % update the Li thickness at the SE/AgC at time it
    LiCC(it+1)       = 10*dt*Li_vmol/FF*XX(Np+2) + LiCC(it);           % update the Li thickness at the CC/AgC at time it    
 
end

%% Data analysis

%% Plot
% 1. Plot the current percentage at the interface and at the particles as a function of time
ttp = linspace(dt,TotT,Nt-1)'/3600;
ifg = ifg + 1;
figure(ifg)
plot(ttp,depSE(2:Nt)/i_ave,'--r',ttp,depCC(2:Nt)/i_ave,'--b');
for ii = 1 : Np
    hold on
    plot(ttp,pi*(dia(ii)/Hx)^2*depI(ii,2:Nt)/i_ave);
end
hold off
xlabel('Charging time (Hours)');
ylabel('Current density rate');
legend('SE/AgC','CC/AgC','Particle1','Particle2');
    
% 2. Plot the voltage as a function of time
ifg = ifg + 1;
figure(ifg)
plot(ttp,-V0(2:Nt)/FF*10^3);
xlabel('Charging time (Hours)');
ylabel('Voltage (mV)');

% 3. Plot Li content in the Ag particles.
ifg = ifg + 1;
figure(ifg)
for ii = 1 : Np
    plot(ttp,yc(ii,2:Nt)*yp_end);
    hold on
end
hold off
xlabel('Charging time (Hours)');
ylabel('y in AgLi_y');
legend('Particle 1','Particle 2');

% Plot Li thickness at the interfaces
ifg = ifg + 1;
figure(ifg)
plot(ttp,LiSE(2:Nt),'--r',ttp,LiCC(2:Nt),'--b');
xlabel('Charging time (Hours)');
ylabel('Li thickness (um)');
legend('SE/AgC','CC/AgC');

%% Define material functions
% The Li chemical potential of AgLi_y vs. Li content (y), in unit J/mol
% chemP_flag=1 means thermo-voltage, chemP_flag=2 means overshoot voltage
function chemPot = chemPot_AgLi(chemP_AgLi_db,chemP_flag,x)
   chemP_db = chemP_AgLi_db(:,chemP_flag+1);
   chemPot = interp1(chemP_AgLi_db(:,1),chemP_db,x);
end

% The mass density of AgLi_y vs. Li content (y), in unit g/cm^3
function massDen = massDensity_AgLi(x)
   massDen = 9.35*exp(-x/2) + 0.57;
end