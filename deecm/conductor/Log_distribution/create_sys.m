%% This script is designed for generating models for distribution of both LPS and NCM particles under different mass ratios of LPS to NCM
%%% Function: create_sys: define the model system
%%% By: Howard Tu, 06/10/2023
function [sys,lps] = create_sys(lps)
%% 1. Define System Parameters in sys structure
sys         = struct;                                         % For massratio study, box size and loading are fixed
sys.load    = 7*10.0;                                                               % The loading is 10mg/cm^2
sys.lx      = 100;                                                                  % The box size in x direction, in unit um 
sys.ly      = 100;                                                                  % The box size in y direction, in unit um
sys.lz      = 150;                                                                  % The box size in z direction, in unit um
sys.prosty  = 1 - sys.load/(sys.lz*lps.den);                                        % The initial porosity
sys.vol     = pi()*sys.lx * sys.ly * sys.lz/4;                                      % The volume of the system
lps.vol    = (1 - sys.prosty) * sys.vol;                                             % The Designed volume of small LPS particles

%% 2. Compute actual particle number and volume
lps = cmp_nmbr(lps,lps.vol);  % For lps particles

%% 3. Plot particle size distribution when debug mode
% Plot theoretical particle size distribution    
% figure(1)
% plot(ncm.dia, 100*ncm.dVol,'-r',ncm.dia, 100*ncm.pVol,'*r');
% hold on
% plot(lps.dia, 100*lps.dVol,'-k',lps.dia, 100*lps.pVol,'*k');
% hold off
% title('Particle distribution')
% xlabel('Diameter (\mum)');
% ylabel('Probability');

%% 4. Other parameters
sys.mass    = lps.rVol*lps.den;                % The total mass of the model, in unit pico-gram
