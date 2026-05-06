%% This script is designed for generating models for distribution of both LPS and NCM particles under different mass ratios of LPS to NCM
%%% Function: create_sys: define the model system
%%% By: Howard Tu, 06/10/2023
function [sys,ncm,lps] = create_sys(ncm, lps,ii,ll)
%% 1. Define System Parameters in sys structure
sys         = struct;                                         % For massratio study, box size and loading are fixed
sys.masRt   = (15-ii)/(80+ii);                                                      % The weight ratio of LPS to NCM
sys.load    = ll*20.0;                                                               % The loading is 10mg/cm^2
sys.lx      = 75;                                                                   % The box size in x direction, in unit um 
sys.ly      = 75;                                                                   % The box size in y direction, in unit um
sys.lz      = 250;                                                                  % The box size in z direction, in unit um
sys.prosty  = 1 - sys.load/(sys.lz*(1+sys.masRt))*(1/ncm.den + sys.masRt/lps.den);  % The initial porosity
sys.vol     = sys.lx * sys.ly * sys.lz;                                             % The volume of the system
sys.volRt   = sys.masRt * ncm.den/lps.den;                                          % The volume ratio of LPS to NCM
lps.vol    = (1 - sys.prosty) * sys.vol * sys.volRt/(1 + sys.volRt);                % The Designed volume of small LPS particles
ncm.vol    = (1 - sys.prosty) * sys.vol * 1/(1 + sys.volRt);                        % The Designed volume of small NCM particles

%% 2. Compute actual particle number and volume
ncm = cmp_nmbr(ncm,ncm.vol);  % For Ag particles 
lps = cmp_nmbr(lps,lps.vol);  % For aC particles

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
sys.mass    = ncm.rVol*ncm.den + lps.rVol*lps.den;                % The total mass of the model, in unit pico-gram
