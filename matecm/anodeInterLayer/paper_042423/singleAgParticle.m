%% This code is developed for the Li-Ag project, for SE-MIEC
% This version considered the Liniearized BV relation, 
% The intraparticle phase separation supressed
% Howard Tu, 07/16/21
clc;clear;
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontName','Arial');
set(0,'defaultAxesFontSize',18);
ifg         = 0;

%% Parameters
% 1. Constants
FF          = 96485.33;              % Faraday constant, in unit s*A/mol
RT          = 8.314*300;             % RT constant, in unit J/mol
Ag_Mmol     = 107.868;               % Molar mass of Silver, in unit g/mol
Ag_Mden     = 9.92;                 % The mass density of Silver, in unit g/cm^3
Ag_Vmol     = Ag_Mmol/Ag_Mden;       % Molar volume of Silver, in unit cm^3/mol
Li_Mmol     = 6.948;                 % Molar mass of Lithium, in unit g/mol
Li_Mden     = 0.534;                 % Mass density of Lithium, in unit g/cm^3
Li_Vmol     = Li_Mmol/Li_Mden;       % Molar volume of Lithium, in unit cm^3/mol
Li_K        = 11;                    % Bulk modulus of Lithium, in unit GPa
% 2. Microstructural parameters
Lx          = 10;                    % Length of the interlayer, in unit um
Hx          = 0.0175;                     % cross-section width, in unit um
BL_por      = 0.5;                   % Inital porosity of the BL layer
Ag_dia      = 0.04;                  % The diameter of Ag particles, in unit um
Np          = length(Ag_dia);        % the number of particles included
% 3. Material properties of the interfaces
i_exc       = 13.0;                   % According to Chiku's paper (Microelectrode studies on charger transfer), unit mA/cm^2
R_Mt        = 1000/(i_exc*FF/RT);    % ASR of charge-transfer resistance at Li-metal/SE interface, in unit Ohm*cm^2
R_SE        = R_Mt/1.0;             % ASR of the charge-transfer resistance at the BL/SS interface, in unit Ohm*cm^2
R_CC        = R_Mt/20.0;             % ASR of the charge-transfer resistance at the BL/CC interface, in unit Ohm*cm^2
R_Pt        = R_Mt/5.0;             % ASR of the charge-transfer resistance at the BL/MetalParticle interface, in unit Ohm*cm^2
% 4. Material properties of the Ag particles
yp_end      = 0.91/(1-0.91);         % The terminal value of y in AgLi_y based on the solubility limit x=0.93
end_Eden    = yp_end*FF/(Ag_Mmol+yp_end*Li_Mmol)*10/36;   % The energy density of AgLi10, in unit mAh/g
end_Mden    = 1.0;                   % Mass density of ending phase AgLi_yend, in unit g/cm^3
end_dia     = Ag_dia .* (Ag_Mden/end_Mden * (1+Li_Mmol/Ag_Mmol*yp_end)).^ (1/3);  % The diameter of end AgLi_y, in unit um
rho         = yp_end*end_Mden/(Ag_Mmol+yp_end*Li_Mmol);    % Interstitial site density of AgLi9, in unit mol/cm^3 = Mass_den / Molar_Mass
% 5. Electrochemical parameters
i_ave       = 0.68;                 % Applied current density, in unit mA/cm^2
 % The chemical potential used:
 % flag=1 means thermo-voltage, flag=2 means overshoot voltage;
 % flag=3 means voltage overshoot to zero from BCC phase directly;
 % flag=4 means voltage overshoot to zero from gamma2 phase
chemP_flag  = 1;                   
% 6. Temperal parameters
Nt          = 200000;                 % The total timesteps
TotT        = 3600;                  % The total discharging time depends on the applied current, in unit s
dt0         = TotT/Nt;               % the time step, in unit s
ttc         = zeros(1,Nt);           % current charging time, in unit s
% 7. Initialization
yc          = zeros(Np,Nt);          % Initialize Li site fraction of each particle, and at the interfaces
yc(:,1)     = 0.01*ones(Np,1);       % The initial Li fraction in Ag particles is 1%
dia         = zeros(Np,Nt);          % The diameter of AgLi_y alloy, in unit um
chem_pot    = zeros(Np,Nt);          % The Li chemical potential of all particles, in unit J/mol
mu_SE       = zeros(1,Nt);           % The Li chemical potential at the SE/AgC interface, in unit J/mol
mu_CC       = zeros(1,Nt);           % The Li chemical potential at the CC/AgC interface, in unit J/mol
depI        = zeros(Np,Nt);          % The deposition current in the particles, in unit mA/cm^2
depSE       = zeros(1,Nt);           % The deposition current at SE/AgC interface, in unit mA/cm^2
depCC       = zeros(1,Nt);           % The deposition current at CC/AgC interface, in unit mA/cm^2
LiSE        = zeros(1,Nt);           % The thickness of deposited Li at the SE/AgC interface, in unit um
LiCC        = zeros(1,Nt);           % The thickness of deposited Li at the CC/AgC interface, in unit um
V0          = zeros(1,Nt);           % The voltage of the cell. In unit J/mol
eqmV        = zeros(Np,Nt);          % The equilibrium voltage of each particles, in unit J/mol
poreLip     = zeros(Np,Nt);           % The volume ratio of Li grow in pores of each AgLi particles
rVolp       = zeros(Np,Nt);           % The volume ratio of each AgLi particles

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
% Entry 1:Np are Li chemical potential in Ag particles; 
% Entry Np+1 is the Li chemical potential at SE/AgC; Entry Np+2 is the Li chemical potential atthe CC/AgC;
% Entry Np+3 is the cell voltage
YY  = [chem_pot(:,1)',mu_SE(1),mu_CC(1),i_ave]';
XX  = zeros(Np+3,1);

% 2. The DD matrix in the equation: yc' = DD*m
DD  = zeros(Np,Np);
cnt = 10*6/(rho*FF);

% 3. Compute and load the Li chemical potential of AgLiy, y from 0 to yp_end, in unit J/mol
% LiAg_freEng_chemPot;
load('chemP_AgLi_db.mat');

%% Numerical solution with Forward Euler method
for it = 1:Nt-1
    if it<Nt/2
        dt = 0.2*dt0;      % Finer timestep in the initial T/10 Hours
    else
        dt = 1.8*dt0;
    end
    
% Step 0: Update the AgLiy particle density and diameters
    AgLi_Mden   = massDensity_AgLi(yc(:,it),yp_end,end_Mden,Ag_Mden,Ag_Mmol,Li_Mmol,Li_Mden);       % The mass density of AgLi_y, in unit g/cm^3
    dia(:,it)   = Ag_dia .* (Ag_Mden ./ AgLi_Mden .* (1+Li_Mmol/Ag_Mmol*yc(:,it)*yp_end)).^ (1/3);  % The diameter of AgLi_y, in unit um
    
    for ip = 1 : Np
       if yc(ip,it) <= 1
            rVolp(ip,it)      = pi/6*dia(ip,it)^3/(Lx*Hx^2);  % The approximated volume ratio of AgLi
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
    chem_pot(:,it)   = chemPot_AgLi(chemP_AgLi_db,chemP_flag,yc(:,it)*yp_end); % The Li chemical potential in Ag particles at current step, in unit J/mol
    YY(1:Np)         = chem_pot(:,it);
    YY(Np+1)         = mu_SE(it);                                       % For Li metal anode, the Li chemical potential is always 0, in unit J/mol
    YY(Np+2)         = mu_CC(it);                                       % For Li metal cathode, the Li chemical potential is always 0, in unit J/mol 
    XX               = MM\YY;                                           % Solve the main equation of all Np+3 unknowns

% Step 2: Compute the accumulate Li content in the Ag particles and at the interfaces
    ytmp           = dt*DD*XX(1:Np) + yc(:,it);                      % update the Li fraction at time it
    SEtmp          = 10*dt*Li_Vmol/FF*XX(Np+1) + LiSE(it);           % update the Li thickness at the SE/AgC at time it
    CCtmp          = 10*dt*Li_Vmol/FF*XX(Np+2) + LiCC(it);           % update the Li thickness at the CC/AgC at time it

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
    LiSE(it+1)       = 10*dt*Li_Vmol/FF*XX(Np+1) + LiSE(it);           % update the Li thickness at the SE/AgC at time it
    LiCC(it+1)       = 10*dt*Li_Vmol/FF*XX(Np+2) + LiCC(it);           % update the Li thickness at the CC/AgC at time it    
 
    ttc(it+1) = ttc(it) + dt;
    
    if ~mod(it,5000)
        it
    end
        
end

%% Data analysis
rVol        = sum(rVolp,1);                     % The volume ratio of AgLi alloy during charging
poreLi      = sum(poreLip,1);                   % The volume ratio of Li in pores in the BL
maxPLi      = BL_por - rVol(end-1);             % The maximal allowed Li volume ratio in the BL
idX         = poreLi>maxPLi;                    % The time when BL is full
extrLi      = zeros(1,Nt);                      % The Li thickness that is extruded out from the BL due to fully dense BL
extrLi(idX) = (poreLi(idX) - maxPLi)*Lx;
poreLi(idX) = maxPLi;
tmpT        = linspace(0,150,1000);              % growth time, in unit second
prs         = FF/Li_Vmol*depI(1,Nt-1)*R_Mt*(1 - exp(-Li_K*Li_Vmol^2*pi*end_dia(1)^2/(R_Mt*FF^2*Hx^2*Lx*BL_por)*tmpT*10^7))/1000; % Pressure, in unit MPa

%% Plot
ttp = linspace(dt,TotT,Nt-1)'/3600;

% 1. Plot the voltage as a function of time
ifg = ifg + 1;
figure(ifg)
plot(ttp*60,-V0(2:Nt)/FF*10^3,'--k');
for ii = 1 : Np
    hold on
    plot(ttp*60,-eqmV(ii,2:Nt)/FF*10^3,'-b');
end
hold off
xlabel('Charging time (Minutes)');
ylabel('Voltage (mV)');
axis([0,1,-100,400]);
legend('AgLi local voltage','AgLi thermo voltage');

% 2. Plot the current percentage at the interface and at the particles as a function of time
ifg = ifg + 1;
figure(ifg)
% yyaxis right
% plot(ttp*60,-V0(2:Nt)/FF*10^3,'--k');
% axis([0,5,-100,400]);
% hold on
% yyaxis left
plot(ttp*60,depSE(2:Nt)/i_ave,'-g',ttp*60,depCC(2:Nt)/i_ave,'-r');
for ii = 1 : Np
    hold on
    plot(ttp*60,depI(ii,2:Nt)/i_ave,'-b');
%     plot(ttp,pi*dia(ii,2:Nt).^2 .* depI(ii,2:Nt)/(i_ave*Hx^2));
end
hold off
xlabel('Charging time (Minutes)');
ylabel('Current density ratio');
axis([0,5,-0.1,1.1]);
legend('SE/BL','CC/BL','AgParticle');

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

% Plot the volume ratio
ifg = ifg + 1;
figure(ifg)
plot(ttp,rVol(1:Nt-1)*100,ttp,poreLi(1:Nt-1)*100);
xlabel('Charging time (Hours)');
ylabel('volume ratio (%)');
legend('AgLi','poreLi');

% Plot Li thickness at the interfaces
ifg = ifg + 1;
figure(ifg)
plot(ttp*3600,LiSE(2:Nt),'--k',ttp*3600,LiCC(2:Nt)+extrLi(2:Nt),'r');
xlabel('Charging time (seconds)');
ylabel('Li thickness (um)');
legend('SE/BL','CC/BL');

%% Define material functions
% The Li chemical potential of AgLi_y vs. Li content (y), in unit J/mol
% chemP_flag=1 means thermo-voltage, chemP_flag=2 means overshoot voltage
function chemPot = chemPot_AgLi(chemP_AgLi_db,chemP_flag,x)
   chemP_db = chemP_AgLi_db(:,chemP_flag+1);
   chemPot = interp1(chemP_AgLi_db(:,1),chemP_db,x);
end

% The mass density of AgLi_y vs. Li content (y), in unit g/cm^3
function massDen = massDensity_AgLi(x,yp_end,end_Mden,Ag_Mden,Ag_Mmol,Li_Mmol,Li_Mden)
   massDen = zeros(length(x),1);
   for ip = 1 : length(x)
     if x(ip) <= 1
         massDen(ip) = (Ag_Mden-end_Mden)*exp(-0.5*x(ip)*yp_end) + end_Mden;
     else
         massDen(ip) = (Ag_Mmol+x(ip)*yp_end*Li_Mmol)/((Ag_Mmol+yp_end*Li_Mmol)/end_Mden + (x(ip)-1)*yp_end*Li_Mmol/Li_Mden);
     end
   end
end