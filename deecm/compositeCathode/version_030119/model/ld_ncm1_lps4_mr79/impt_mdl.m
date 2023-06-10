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
ncm.min     = ncm.ave / ncm.ave(1) * 3.2;   % The minimum diameter of NMC particles, reference to 5um-NMC
ncm.max     = ncm.ave / ncm.ave(1) * 8.5;   % The maximum diameter of NMC particles, reference to 5um-NMC
ncm.sgm     = ncm.ave / ncm.ave(1) * 1.0;   % The deviation of NMC particles, reference to 5um-NMC

% Parameters for each NMC size
for itmp = 1 : par.nNcm
    stmp = strcat('ncm.s',num2str(itmp),'.min = ncm.min(',num2str(itmp),');');
    eval(stmp);                                % The minimum diameter of current NMC size
    stmp = strcat('ncm.s',num2str(itmp),'.max = ncm.max(',num2str(itmp),');');
    eval(stmp);                                % The maximum diameter of current NMC size
    stmp = strcat('ncm.s',num2str(itmp),'.ave = ncm.ave(',num2str(itmp),');');
    eval(stmp);                                % The average diameter of current NMC size
    stmp = strcat('ncm.s',num2str(itmp),'.sgm = ncm.sgm(',num2str(itmp),');');
    eval(stmp);
end

%% 2.Define properties of different sizes LPS particles in lps structure
lps         = struct;
lps.nD      = 20;            % The number of sampling diameters of LPS particles
lps.den     = 1.87;          % The density of LPS material, in unit picogram/micron^3
lps.typ     = 2;             % The type of LPS particles
lps.ave     = [0.4,0.8,1.2,1.5,3.0,5,8,12,15];   % The designed average diameter of all LPS, in unit um
lps.min     = lps.ave / lps.ave(4) * 0.95;       % The minimum diameter of LPS particles, reference to 1.5um-LPS
lps.max     = lps.ave / lps.ave(4) * 3.1;        % The maximum diameter of LPS particles, reference to 1.5um-LPS
lps.sgm     = lps.ave / lps.ave(4) * 0.5;        % The deviation of LPS particles, reference to 1.5um-LPS

% Parameters for each LPS size
for itmp = 1 : par.nLps
    stmp = strcat('lps.s',num2str(itmp),'.min = lps.min(',num2str(itmp),');');
    eval(stmp);                                % The minimum diameter of current LPS size
    stmp = strcat('lps.s',num2str(itmp),'.max = lps.max(',num2str(itmp),');');
    eval(stmp);                                % The maximum diameter of current LPS size
    stmp = strcat('lps.s',num2str(itmp),'.ave = lps.ave(',num2str(itmp),');');
    eval(stmp);                                % The average diameter of current LPS size
    stmp = strcat('lps.s',num2str(itmp),'.sgm = lps.sgm(',num2str(itmp),');');
    eval(stmp);                                % The deviation of current LPS size
end

%% 3.Define ghost particles at the top of the system to apply pressure
gst         = struct;
gst.den     = 10;            % The density of ghost material
gst.typ     = 3;

%% 4. Compute theoretical particle probabilities of all particles
% For NMC particles
for itmp = 1 : par.nNcm
    tmp = strcat('ncm.s',num2str(itmp),'  = cmp_prob(ncm.s',num2str(itmp),',ncm.nD,1);');
    eval(tmp);
end
% For LPS particles
for itmp = 1 : par.nLps
    tmp = strcat('lps.s',num2str(itmp),'  = cmp_prob(lps.s',num2str(itmp),',lps.nD,1);');
    eval(tmp);
end
