%% This script is designed for generating models for distribution of both Ag and active carbon particles
%%% Function: mat_Inpt: define geometry and different types of material and their properties
%%% By: Howard Tu, 05/28/2023
function [sys,ag,ac] = create_sys(ii, ag, ac)
%% 1. Define System Parameters in sys structure
sys         = struct;                                         % For massratio study, box size and loading are fixed
sys.masRt   = 1.0/ii;                                         % The weight ratio of Ag to C
sys.lx      = 2;                                              % The box size in x direction, in unit um 
sys.ly      = 2;                                              % The box size in y direction, in unit um
sys.lzf     = 10;                                             % The desired box size in z direction after packing, in unit um
% The desired porosity of the interlayer after packing, can be obtained from characterization; 
% however, the theoretical smallest value is 0.26 for HCP pack, and 0.365 for random pack
sys.prf     = 0.4;                                            
sys.lz      = sys.lzf*3;                                      % The box size in z direction before packing.
sys.pr      = 1-(1-sys.prf)*sys.lzf/sys.lz;                   % The porosity of the interlayer before packing
sys.vol     = sys.lx * sys.ly * sys.lz;                       % The volume of the system
sys.volRt   = sys.masRt * ac.den/ag.den;                      % The volume ratio of Ag to C
ag.vol      = (1 - sys.pr) * sys.vol * sys.volRt/(1 + sys.volRt);    % The designed volume of all Ag particles
ac.vol      = ag.vol/sys.volRt;                               % The designed volume of all active carbon particles

%% 2. Compute actual particle number and volume
% For Ag particles
ag.nTot     = round(6/pi * ag.vol/ag.dia^3);                  % The number of Ag particles
ag.rVol     = pi/6 * ag.dia^3 * ag.nTot;                      % The real volume of all Ag particles

% For active carbon particles
ac.nTot     = round(6/pi * ac.vol/ac.dia^3);                  % The number of C particles
ac.rVol     = pi/6 * ac.dia^3 * ac.nTot;                      % The real volume of all C particles

% Other parameters
sys.mass    = ag.rVol*ag.den + ac.rVol*ac.den;                % The total mass of the model, in unit pico-gram
