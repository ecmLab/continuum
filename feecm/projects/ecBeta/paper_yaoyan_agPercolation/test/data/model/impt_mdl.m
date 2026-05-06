%% This script is designed for generating models for distribution of both LPS and NCM particles
%%% Function: impt_mat: define geometry and material properties of all particles
%%% By: Howard Tu
function [ncm,lps,gst] = impt_mdl(par)
%% 1.Define properties of NCM particles in ncm structure
ncm         = struct;
ncm.nD      = 20;            % The number of sampling diameters of NCM particles
ncm.den     = 4.85;          % The density of NCM material, in unit picogram/micron^3
ncm.typ     = 1;             % The type of NCM particles
ncm.ave     = [5.0, 10, 12];                % The designed average diameter of all NMC, in unit um

% Parameters for each NMC size
for itmp = 1 : par.nNcm
    stmp = strcat('ncm.s',num2str(itmp),'.dia = ncm.ave(',num2str(itmp),');');
    eval(stmp);                                % The average diameter of current NMC size     
end

%% 2.Define properties of different sizes LPS particles in lps structure
lps         = struct;
lps.nD      = 20;            % The number of sampling diameters of LPS particles
lps.den     = 1.87;          % The density of LPS material, in unit picogram/micron^3
lps.typ     = 2;             % The type of LPS particles
lps.ave     = [0.4,1.0,2.0,3.0,4.0,5.0,6.0,7.0];   % The designed average diameter of all LPS, in unit um

% Parameters for each LPS size
for itmp = 1 : par.nLps
    stmp = strcat('lps.s',num2str(itmp),'.dia = lps.ave(',num2str(itmp),');');
    eval(stmp);                                % The average diameter of current LPS size
end

%% 3.Define ghost particles at the top of the system to apply pressure
gst         = struct;
gst.den     = 10;            % The density of ghost material
gst.typ     = 3;

