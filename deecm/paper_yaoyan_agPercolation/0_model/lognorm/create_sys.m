%% This script is designed for generating models for distribution of both Ag and active carbon particles
%%% Function: create_sys: define the model system
%%% By: Howard Tu, 06/10/2023
function [sys,ag,ac] = create_sys(ag, ac)
%% 1. Define System Parameters in sys structure
sys         = struct;                                         % For massratio study, box size and loading are fixed
sys.masRt   = 1.0/3.0;                                         % The weight ratio of Ag to C
sys.lx      = 1.0;                                              % The box size in x direction, in unit um; 20 is for plot; 2 is for calculation
sys.ly      = 1.0;                                              % The box size in y direction, in unit um; 20 is for plot; 2 is for calculation
sys.lzf     = 10;                                             % The desired box size in z direction after packing, in unit um
% The desired porosity of the interlayer after packing, can be obtained from characterization; 
% however, the theoretical smallest value is 0.26 for HCP pack, and 0.365 for random pack
sys.prf     = 0.4;                                            
sys.lz      = sys.lzf*2;                                      % The box size in z direction before packing.
sys.pr      = 1-(1-sys.prf)*sys.lzf/sys.lz;                   % The porosity of the interlayer before packing
sys.vol     = sys.lx * sys.ly * sys.lz;                       % The volume of the system
sys.volRt   = sys.masRt * ac.den/ag.den;                      % The volume ratio of Ag to C
ag.vol      = (1 - sys.pr) * sys.vol * sys.volRt/(1 + sys.volRt);    % The designed volume of all Ag particles
ac.vol      = ag.vol/sys.volRt;                               % The designed volume of all active carbon particles

%% 2. Compute actual particle number and volume
ag = cmp_nmbr(ag,ag.vol);  % For Ag particles 
ac = cmp_nmbr(ac,ac.vol);  % For aC particles

%% 3. Plot particle size distribution when debug mode
% Plot theoretical particle size distribution    
figure(1)
plot(ag.dia, 100*ag.dVol,'-r',ag.dia, 100*ag.pVol,'*r');
hold on
plot(ac.dia, 100*ac.dVol,'-k',ac.dia, 100*ac.pVol,'*k');
hold off
title('Particle distribution')
xlabel('Diameter (\mum)');
ylabel('Probability');

%% 4. Other parameters
sys.mass    = ag.rVol*ag.den + ac.rVol*ac.den;                % The total mass of the model, in unit pico-gram
