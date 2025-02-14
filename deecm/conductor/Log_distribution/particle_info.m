%% This script is designed for generating models for distribution of both NMC and C particles
%%% Function: impt_mat: define geometry and material properties of all particles
%%% By: Howard Tu, 06/10/2023
function [lps,gst] = particle_info(d)
%% 1.Define properties of active LPS particles in LPS structure
lps         = struct;
lps.nD      = 29;            % The number of sampling diameters of active LPS particles
lps.den     = 1.87;           % The density of active LPS material, in unit picogram/micron^3, equivalent to g/cm^3
lps.typ     = 1;             % Label active LPS particles as Type 2 particles
lps.ave     = 3;          % The designed diameter of active LPS particles, in unit um; 0.2 is for plot; 0.02 is for calculation
lps.min     = 1;         % The minimum diameter of LPS particles
lps.max     = 5;          % The maximum diameter of LPS particles
lps.sgm     = d;          % The deviation of LPS particles

%% 3.Define ghost particles at the top of the system to apply pressure
gst         = struct;
gst.den     = 10;            % The density of ghost material
gst.typ     = 2;

%% 4. Compute theoretical particle probabilities of all particles
lps = cmp_prob(lps,lps.nD,1);  % For LPS particles