%% This script is designed for generating models for distribution of both Ag and C particles
%%% Function: impt_mat: define geometry and material properties of all particles
%%% By: Howard Tu, 06/10/2023
function [ag,ac,gst] = particle_info(ii)
%% 1.Define properties of Ag particles in ag structure
ag          = struct;
ag.nD       = 20;            % The number of sampling diameters of Ag particles
ag.den      = 10.49;         % The density of Ag material, in unit picogram/micron^3, equivalent to g/cm^3
ag.typ      = 1;             % Label Ag particles as Type 1 particles
if ii == 1           % small size + fine mix: Ag_ave = 40nm, narrow lognorm distribution
    ag.ave     = 0.04;       % The designed average diameter of all Ag particles, in unit um
    ag.min     = 0.035;      % The minimum diameter of Ag particles
    ag.max     = 0.045;       % The maximum diameter of Ag particles
    ag.sgm     = 0.002;       % The deviation of Ag particles
elseif ii == 2   % small size + hand mix: Ag_ave = 40nm, wide lognorm distribution
    ag.ave     = 0.04;       % The designed average diameter of all Ag particles, in unit um
    ag.min     = 0.03;      % The minimum diameter of Ag particles
    ag.max     = 0.08;       % The maximum diameter of Ag particles
    ag.sgm     = 0.01;       % The deviation of Ag particles
elseif ii == 3   % large size + fine mix: Ag_ave = 80nm, narrow lognorm distribution
    ag.ave     = 0.08;       % The designed average diameter of all Ag particles, in unit um
    ag.min     = 0.075;      % The minimum diameter of Ag particles
    ag.max     = 0.085;       % The maximum diameter of Ag particles
    ag.sgm     = 0.002;       % The deviation of Ag particles
elseif ii == 4   % large size + hand mix: Ag_ave = 80nm, wide lognorm distribution
    ag.ave     = 0.08;       % The designed average diameter of all Ag particles, in unit um
    ag.min     = 0.065;      % The minimum diameter of Ag particles
    ag.max     = 0.16;       % The maximum diameter of Ag particles
    ag.sgm     = 0.02;       % The deviation of Ag particles
end

%% 2.Define properties of active carbon particles in ac structure
ac         = struct;
ac.nD      = 20;            % The number of sampling diameters of active carbon particles
ac.den     = 2.1;           % The density of active carbon material, in unit picogram/micron^3, equivalent to g/cm^3
ac.typ     = 2;             % Label active carbon particles as Type 2 particles
ac.ave     = 0.02;          % The designed diameter of active carbon particles, in unit um; 0.2 is for plot; 0.02 is for calculation
ac.min     = 0.015;         % The minimum diameter of ac particles
ac.max     = 0.04;          % The maximum diameter of ac particles
ac.sgm     = 0.005;          % The deviation of ac particles

%% 3.Define ghost particles at the top of the system to apply pressure
gst         = struct;
gst.den     = 10;            % The density of ghost material
gst.typ     = 3;

%% 4. Compute theoretical particle probabilities of all particles
ag = cmp_prob(ag,ag.nD,1);  % For Ag particles 
ac = cmp_prob(ac,ac.nD,1);  % For aC particles