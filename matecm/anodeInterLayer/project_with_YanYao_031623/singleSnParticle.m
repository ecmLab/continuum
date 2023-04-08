%% This code is developed for the Li-Sn project, for SE-MIEC
% This version considered the Liniearized BV relation, 
% The intraparticle phase separation supressed
% Howard Tu, 07/26/21
clc;clear;
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontName','Arial');
set(0,'defaultAxesFontSize',18);
ifg         = 0;

%% Parameters
% 1. Constants
FF          = 96485.33;              % Faraday constant, in unit s*A/mol
RT          = 8.314*300;             % RT constant, in unit J/mol
Sn_Mmol     = 118.71;               % Molar mass of Silver, in unit g/mol
Sn_Mden     = 7.31;                 % The mass density of Silver, in unit g/cm^3
Sn_Vmol     = Sn_Mmol/Sn_Mden;       % Molar volume of Silver, in unit cm^3/mol
Li_Mmol     = 6.948;                 % Molar mass of Lithium, in unit g/mol
Li_Mden     = 0.534;                 % Mass density of Lithium, in unit g/cm^3
Li_Vmol     = Li_Mmol/Li_Mden;       % Molar volume of Lithium, in unit cm^3/mol
Li_K        = 11;                    % Bulk modulus of Lithium, in unit GPa
% 2. Microstructural parameters
Lx          = 10;                    % Length of the interlayer, in unit um
Hx          = 2;                     % cross-section width, in unit um
BL_por      = 0.5;                   % Inital porosity of the BL layer
Sn_dia      = 0.95;                  % The diameter of Sn particles, in unit um
Np          = length(Sn_dia);        % the number of particles included
% 3. Material properties of the interfaces
i_exc       = 0.3;                   % According to Chiku's paper (Microelectrode studies on charger transfer), unit mA/cm^2
R_Mt        = 1000/(i_exc*FF/RT);    % ASR of charge-transfer resistance at Li-metal/SE interface, in unit Ohm*cm^2
R_SE        = R_Mt/1.0;             % ASR of the charge-transfer resistance at the BL/SS interface, in unit Ohm*cm^2
R_CC        = R_Mt/1.0;             % ASR of the charge-transfer resistance at the BL/CC interface, in unit Ohm*cm^2
R_Pt        = R_Mt/1.0;             % ASR of the charge-transfer resistance at the BL/MetalParticle interface, in unit Ohm*cm^2
% 4. Material properties of the Sn particles
yp_end      = 4.4;                  % The terminal value of y in SnLi_y is SnLi7
end_Eden    = yp_end*FF/(Sn_Mmol+yp_end*Li_Mmol)*10/36;   % The energy density of SnLi7, in unit mAh/g
end_Mden    = 2.583;                   % Mass density of ending phase SnLi_yend, in unit g/cm^3
end_dia     = Sn_dia .* (Sn_Mden/end_Mden * (1+Li_Mmol/Sn_Mmol*yp_end)).^ (1/3);  % The diameter of end SnLi_y, in unit um
rho         = yp_end*end_Mden/(Sn_Mmol+yp_end*Li_Mmol);    % Interstitial site density of SnLi9, in unit mol/cm^3 = Mass_den / Molar_Mass
% 5. Electrochemical parameters
i_ave       = 0.68;                 % Applied current density, in unit mA/cm^2
% 6. Temperal parameters
Nt          = 100000;                 % The total timesteps
TotT        = 3600;                  % The total discharging time depends on the applied current, in unit s
dt          = TotT/Nt;               % the time step, in unit s
% 7. Initialization
yc          = zeros(Np,Nt);          % Initialize Li site fraction of each particle, and at the interfaces
yc(:,1)     = 0.01*ones(Np,1);       % The initial Li fraction in Sn particles is 2%
dia         = zeros(Np,Nt);          % The diameter of SnLi_y alloy, in unit um
chem_pot    = zeros(Np,Nt);          % The Li chemical potential of all particles, in unit J/mol
mu_SE       = zeros(1,Nt);           % The Li chemical potential at the SE/SnC interface, in unit J/mol
mu_CC       = zeros(1,Nt);           % The Li chemical potential at the CC/SnC interface, in unit J/mol
depI        = zeros(Np,Nt);          % The deposition current in the particles, in unit mA/cm^2
depSE       = zeros(1,Nt);           % The deposition current at SE/SnC interface, in unit mA/cm^2
depCC       = zeros(1,Nt);           % The deposition current at CC/SnC interface, in unit mA/cm^2
LiSE        = zeros(1,Nt);           % The thickness of deposited Li at the SE/SnC interface, in unit um
LiCC        = zeros(1,Nt);           % The thickness of deposited Li at the CC/SnC interface, in unit um
V0          = zeros(1,Nt);           % The voltage of the cell. In unit J/mol
eqmV        = zeros(Np,Nt);          % The equilibrium voltage of each particles, in unit J/mol
poreLip     = zeros(Np,Nt);           % The volume ratio of Li grow in pores of each SnLi particles
rVolp       = zeros(Np,Nt);           % The volume ratio of each SnLi particles

%% Compute all the matrix constant for solving equations
% 1. Components of MM in the matrix equation: MM*X = Y
M11            = -FF*R_Pt/10^3*eye(Np+2);  % In unit J/mol/(mA/cm^2)
M11(Np+1,Np+1) = -FF*R_SE/10^3;
M11(Np+2,Np+2) = -FF*R_CC/10^3;
M12            = ones(Np+2,1);             % in unit 1
M21            = ones(1,Np+2);             % In unit 1
M22            = 0;
% Assemble the MM matrix, size of MM: (Np+3) x (Np+3)
MM  = [M11, M12; M21, M22];
% Defind Y in the matrix equation: MM*X = Y. Size of Y vector: (Np+3) x 1
% Entry 1:Np are Li chemical potential in Sn particles; 
% Entry Np+1 is the Li chemical potential at SE/SnC; Entry Np+2 is the Li chemical potential atthe CC/SnC;
% Entry Np+3 is the cell voltage
YY  = [chem_pot(:,1)',mu_SE(1),mu_CC(1),i_ave]';
XX  = zeros(Np+3,1);

% 2. The DD matrix in the equation: yc' = DD*m
DD  = zeros(Np,Np);
cnt = 10*6/(rho*FF);

% 3. Compute and load the Li chemical potential of SnLiy, y from 0 to yp_end, in unit J/mol
load('chemP_SnLi_db.mat');

%% Numerical solution with Forward Euler method
for it = 1:Nt-1

% Step 0: Update the SnLiy particle density and diameters
    SnLi_Mden   = massDensity_SnLi(yc(:,it),yp_end,end_Mden,Sn_Mden);       % The mass density of SnLi_y, in unit g/cm^3
    dia(:,it)   = Sn_dia .* (Sn_Mden ./ SnLi_Mden .* (1+Li_Mmol/Sn_Mmol*yc(:,it)*yp_end)).^ (1/3);  % The diameter of SnLi_y, in unit um
    
    for ip = 1 : Np
       if yc(ip,it) <= 1
            rVolp(ip,it)      = pi/6*dia(ip,it)^3/(Lx*Hx^2);  % The approximated volume ratio of SnLi
       else
            rVolp(ip,it)      = pi/6*end_dia(ip)^3/(Lx*Hx^2);
            poreLip(ip,it)    = pi/6*(dia(ip,it)^3-end_dia(ip)^3)/(Lx*Hx^2);  % The approximated volume ratio of Li in pores
       end
    end
    
    for ip = 1 : Np
        M21(ip)   = pi/Hx^2 * dia(ip,it).^2;
        DD(ip,ip) = cnt/dia(ip,it);       % In unit 1/(mA/cm^2) * (1/s)
    end

% Step 1: Try the first calculation
    chem_pot(:,it)   = chemPot_SnLi(chemP_SnLi_db,yc(:,it)*yp_end); % The Li chemical potential in Sn particles at current step, in unit J/mol
    YY(1:Np)         = chem_pot(:,it);
    YY(Np+1)         = mu_SE(it);                                       % For Li metal anode, the Li chemical potential is always 0, in unit J/mol
    YY(Np+2)         = mu_CC(it);                                       % For Li metal cathode, the Li chemical potential is always 0, in unit J/mol 
    XX               = MM\YY;                                           % Solve the main equation of all Np+3 unknowns

% Step 2: Compute the accumulate Li content in the Sn particles and at the interfaces
    ytmp           = dt*DD*XX(1:Np) + yc(:,it);                      % update the Li fraction at time it
    SEtmp          = 10*dt*Li_Vmol/FF*XX(Np+1) + LiSE(it);           % update the Li thickness at the SE/SnC at time it
    CCtmp          = 10*dt*Li_Vmol/FF*XX(Np+2) + LiCC(it);           % update the Li thickness at the CC/SnC at time it

% Step 3: Impose constraints: 
% (1) Li content in particles cannot be less than zero; 
% (2) Li thickness at interfaces cannot be less than zero
    [tmp,yc0]      = find(ytmp<0);
    [tmp,yct]      = find([SEtmp,CCtmp]<0);
    if (length(yc0)+length(yct)) > 0
       MMt         = MM;
       YYt         = YY;
       for iy = 1:length(yc0)
          MMt(yc0(iy),yc0(iy)) = 1;
          MMt(yc0(iy),Np+3)    = 0;
          YYt(yc0(iy))         = 0;
       end

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
    eqmV(:,it+1)     = YY(1:Np);
    yc(:,it+1)       = dt*DD*depI(:,it+1) + yc(:,it);
    LiSE(it+1)       = 10*dt*Li_Vmol/FF*XX(Np+1) + LiSE(it);           % update the Li thickness at the SE/SnC at time it
    LiCC(it+1)       = 10*dt*Li_Vmol/FF*XX(Np+2) + LiCC(it);           % update the Li thickness at the CC/SnC at time it    
 
end

%% Data analysis
rVol        = sum(rVolp,1);                     % The volume ratio of SnLi alloy during charging
poreLi      = sum(poreLip,1);                   % The volume ratio of Li in pores in the BL
maxPLi      = BL_por - rVol(end-1);             % The maximal allowed Li volume ratio in the BL
idX         = poreLi>maxPLi;                    % The time when BL is full
extrLi      = zeros(1,Nt);                      % The Li thickness that is extruded out from the BL due to fully dense BL
extrLi(idX) = (poreLi(idX) - maxPLi)*Lx;
poreLi(idX) = maxPLi;
tmpT        = linspace(0,150,1000);              % growth time, in unit second
prs         = FF/Li_Vmol*depI(1,Nt-1)*R_Mt*(1 - exp(-Li_K*Li_Vmol^2*pi*end_dia(1)^2/(R_Mt*FF^2*Hx^2*Lx*BL_por)*tmpT*10^7))/1000; % Pressure, in unit MPa

%% Plot
% 1. Plot the current percentage at the interface and at the particles as a function of time
ttp = linspace(dt,TotT,Nt-1)'/3600;
ifg = ifg + 1;
figure(ifg)
plot(ttp,depSE(2:Nt)/i_ave,'--r',ttp,depCC(2:Nt)/i_ave,'--b');
for ii = 1 : Np
    hold on
    plot(ttp,depI(ii,2:Nt)/i_ave);
%     plot(ttp,pi*dia(ii,2:Nt).^2 .* depI(ii,2:Nt)/(i_ave*Hx^2));
end
hold off
xlabel('Charging time (Hours)');
ylabel('Current density rate');
legend('SE/BL','CC/BL','SnParticle');
    
% 2. Plot the voltage as a function of time
ifg = ifg + 1;
figure(ifg)
plot(ttp,-V0(2:Nt)/FF*10^3);
for ii = 1 : Np
    hold on
    plot(ttp,-eqmV(ii,2:Nt)/FF*10^3);
end
plot(chemP_SnLi_db(:,1)/70,-chemP_SnLi_db(:,2)/FF*10^3);
hold off
xlabel('Charging time (Hours)');
ylabel('Voltage (mV)');
legend('Cell voltage','Equilibrium voltage');

% 3. Plot Li content in the Sn particles.
ifg = ifg + 1;
figure(ifg)
for ii = 1 : Np
    plot(ttp,yc(ii,2:Nt)*yp_end);
    hold on
end
hold off
xlabel('Charging time (Hours)');
ylabel('y in SnLi_y');

% Plot the volume ratio
ifg = ifg + 1;
figure(ifg)
plot(ttp,rVol(1:Nt-1)*100,ttp,poreLi(1:Nt-1)*100);
xlabel('Charging time (Hours)');
ylabel('volume ratio (%)');
legend('SnLi','poreLi');

% Plot Li thickness at the interfaces
ifg = ifg + 1;
figure(ifg)
plot(ttp,LiSE(2:Nt),'--r',ttp,LiCC(2:Nt)+extrLi(2:Nt),'--b');
xlabel('Charging time (Hours)');
ylabel('Li thickness (um)');
legend('SE/BL','CC/BL');

% Plot the density change
ifg = ifg + 1;
figure(ifg)
xx   =linspace(0,2,100);
yy   = massDensity_SnLi(xx,yp_end,end_Mden,Sn_Mden);
Lixx =[0.3333    0.4000    1.0000    2.3333    2.5000    2.6000    3.0000    3.5000    4.2500    yp_end];
royy =[5.9950    5.9670    4.9880    3.6540    3.5560    3.4840    3.2600    2.9750    2.5900    end_Mden];
plot(xx*yp_end,yy, Lixx, royy, 'o');
xlabel('y in SnLiy');
ylabel('mass density (g/cm^3)');

% Plot the internal pressure
ifg = ifg + 1;
figure(ifg)
plot(tmpT, prs);
xlabel('Time (s)');
ylabel('Pressure (MPa)');

%% Define material functions
% The Li chemical potential of SnLi_y vs. Li content (y), in unit J/mol
% chemP_flag=1 means thermo-voltage, chemP_flag=2 means overshoot voltage
function chemPot = chemPot_SnLi(chemP_SnLi_db,x)
   if x < chemP_SnLi_db(end,1)
        chemPot = interp1(chemP_SnLi_db(:,1),chemP_SnLi_db(:,2),x);
   else
        chemPot = 0;
   end
end

% The mass density of SnLi_y vs. Li content (y), in unit g/cm^3
function massDen = massDensity_SnLi(x,yp_end,end_Mden,Sn_Mden)
   massDen = zeros(length(x),1);
   for ip = 1 : length(x)
       massDen(ip) = (Sn_Mden-end_Mden)*exp(-0.7*x(ip)*yp_end) + end_Mden;
   end
end