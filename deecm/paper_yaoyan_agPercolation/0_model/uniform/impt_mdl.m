%% This script is designed for generating models for distribution of both Ag and C particles
%%% Function: impt_mat: define geometry and material properties of all particles
%%% By: Howard Tu, 05/28/2023

%% 1.Define properties of NCM particles in ncm structure
ag          = struct;
ag.nD       = 20;            % The number of sampling diameters of Ag particles
ag.den      = 10.49;         % The density of Ag material, in unit picogram/micron^3, equivalent to g/cm^3
ag.typ      = 1;             % Label Ag particles as Type 1 particles
ag.dia      = 0.8;           % The designed diameter of all NMC particles, in unit um; 0.8 is for plot; 0.1 is for calculation

%% 2.Define properties of different sizes active carbon particles in lps structure
ac         = struct;
ac.nD      = 20;            % The number of sampling diameters of active carbon particles
ac.den     = 2.1;          % The density of active carbon material, in unit picogram/micron^3, equivalent to g/cm^3
ac.typ     = 2;             % Label active carbon particles as Type 2 particles
ac.dia     = 0.4;           % The designed diameter of active carbon particles, in unit um; 0.4 is for plot; 0.05 is for calculation

%% 3.Define ghost particles at the top of the system to apply pressure
gst         = struct;
gst.den     = 10;            % The density of ghost material
gst.typ     = 3;

