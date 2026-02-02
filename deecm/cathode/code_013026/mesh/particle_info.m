%% This script is designed for generating models for distribution of both NMC and C particles
%%% Function: impt_mat: define geometry and material properties of all particles
%%% By: Howard Tu, 06/10/2023
function [ncm,lps,gst] = particle_info(iw)
%% 1.Define properties of NMC particles in NMC structure
ncm         = struct;
ncm.nD      = 29;            % The number of sampling diameters of NMC particles
ncm.den     = 4.85;          % The density of NMC material, in unit picogram/micron^3, equivalent to g/cm^3
ncm.typ     = 1;             % Label NMC particles as Type 1 particles
ncm.ave     = 5;       % The designed averNMCe diameter of all NMC particles, in unit um
ncm.min     = 3.5;         % The minimum diameter of NMC particles
ncm.max     = 6.5;         % The maximum diameter of NMC particles
ncm.sgm     = 0.1*iw;       % The deviation of NMC particles


%% 2.Define properties of active LPS particles in LPS structure
lps         = struct;
lps.nD      = 29;            % The number of sampling diameters of active LPS particles
lps.den     = 1.87;           % The density of active LPS material, in unit picogram/micron^3, equivalent to g/cm^3
lps.typ     = 2;             % Label active LPS particles as Type 2 particles
lps.ave     = 2.5;          % The designed diameter of active LPS particles, in unit um; 0.2 is for plot; 0.02 is for calculation
lps.min     = 1;         % The minimum diameter of LPS particles
lps.max     = 4;          % The maximum diameter of LPS particles
lps.sgm     = 0.1;          % The deviation of LPS particles

%% 3.Define ghost particles at the top of the system to apply pressure
gst         = struct;
gst.den     = 10;            % The density of ghost material
gst.typ     = 3;

%% 4. Compute theoretical particle probabilities of all particles
ncm = cmp_prob(ncm,ncm.nD,1);  % For NMC particles 
lps = cmp_prob(lps,lps.nD,1);  % For LPS particles