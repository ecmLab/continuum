%% Create System for DRX-Carbon Mixture
%%% Function: create_sys: define the model system for DRX (active) and carbon (conductive) particles
%%% By: Howard Tu, 06/10/2023 (mod. DRX/C terminology)
function [sys,drx,carbon] = create_sys(drx, carbon, ii)
%% 1. Define System Parameters in sys structure
sys         = struct;                                         % For mass-ratio study, box size and loading are fixed
sys.masRt   = (6.5-ii*0.5)/(88.5+ii*0.5);                       % The weight ratio of carbon to DRX, assuming 5 wt% of binder
                                                              % Carbon wt% changes from 6% to 1%
sys.load    = 10.0;                                           % The areal loading for DRX is generally load, < 10mg/cm^2, according to:
                                                              % https://www.sciencedirect.com/science/article/pii/S2542435121005316
sys.lx      = 15;                                             % Box size in x direction [μm]
sys.ly      = 15;                                             % Box size in y direction [μm]
sys.lz      = 30;                                            % Box size in z direction [μm]
sys.prosty  = 1 - sys.load/(sys.lz*(1+sys.masRt))*(1/drx.den + sys.masRt/carbon.den);  % Initial porosity
sys.vol     = sys.lx * sys.ly * sys.lz;                       % Volume of the system
sys.volRt   = sys.masRt * drx.den/carbon.den;                 % Volume ratio of carbon to DRX
carbon.vol  = (1 - sys.prosty) * sys.vol * sys.volRt/(1 + sys.volRt);  % Designed volume of carbon particles
drx.vol     = (1 - sys.prosty) * sys.vol * 1/(1 + sys.volRt);         % Designed volume of DRX particles

%% 2. Compute actual particle number and volume
drx   = cmp_nmbr(drx,drx.vol);      % For DRX particles

% Carbon clusters (SuperP aggregates of 0.1 μm spheres)
sphere_volume = (pi/6) * carbon.dia^3;
cluster_volume = carbon.cluster.cluster_volume;
carbon.nClusters = max(1, round(carbon.vol / cluster_volume));
carbon.nTot      = carbon.nClusters * carbon.cluster.sphere_count;
carbon.rVol      = carbon.nTot * sphere_volume;
carbon.nDia      = carbon.nTot;

%% 3. Plot particle size distribution when debug mode
%Plot theoretical particle size distribution if needed
figure(1)
plot(drx.dia, 100*drx.dVol,'-r',drx.dia, 100*drx.pVol,'*r');
%hold on
%plot(carbon.dia, 100*carbon.dVol,'-k',carbon.dia, 100*carbon.pVol,'*k');
hold off
title('Particle distribution')
xlabel('Diameter (\mum)');
ylabel('Probability');

%% 4. Other parameters
sys.mass    = drx.rVol*drx.den + carbon.rVol*carbon.den;      % Total mass [pg]
