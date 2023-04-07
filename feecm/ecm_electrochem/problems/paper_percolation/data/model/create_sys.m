%% This script is designed for generating models for distribution of both LPS and NCM particles
%%% Function: mat_Inpt: define geometry and different types of material and their properties
%%% By: Howard Tu
function [sys,ncm,lps] = create_sys(ii, par, ncm, lps)
%% 1. Define System Parameters in sys structure
sys             = struct;

if par.iSty     == 1                                      % For box size study, massratio, loading and thickness are fixed
    sys.masRt   = (95-par.ncmMr)/par.ncmMr;               % The weight ratio of LPS to NCM if fixed, 56,80, 90 are studied
% Scale the system equivalently according to the NMC particle size used in the current study
    sys.load    = 10*10.0*ncm.ave(par.iNcm)/ncm.ave(1);   % The loading is 10mg/cm^2 at 5um-NMC, other NMC are scaled accordingly
    sys.lz      = 80*ncm.ave(par.iNcm)/ncm.ave(1);        % The initial thickness is 80um at 5um-NMC, other NMC are scaled accordingly
    sys.lxyMn   = 30*ncm.ave(par.iNcm)/ncm.ave(1);        % The minimal length of intersection is  30um at 5um-NMC, other NMC are scaled accordingly
    sys.lxyMx   = 120*ncm.ave(par.iNcm)/ncm.ave(1);       % The maximal length of intersection is 120um at 5um-NMC, other NMC are scaled accordingly
    lxy         = sys.lxyMn*(sys.lxyMx/sys.lxyMn)^((ii-1)/(par.nMdl-1));                   % Area length increases as a geometric sequence
    sys.lx      = lxy;                                    % The box size in x direction
    sys.ly      = lxy;                                    % The box size in y direction
    sys.prosty  = 1 - sys.load/(sys.lz*(1+sys.masRt))*(1/ncm.den + sys.masRt/lps.den); % The initial porosity
    
elseif par.iSty == 2                                      % For loading study, massratio, initial porosity are fixed
    sys.masRt   = (95-par.ncmMr)/par.ncmMr;               % The weight ratio of LPS to NCM if fixed, 56,80, 90 are studied
    sys.prosty  = 0.8;                                    % The initial porosity
% Scale the system equivalently according to the NMC particle size used in the current study
    sys.lx      = 80*ncm.ave(par.iNcm)/ncm.ave(1);        % The box size in x direction is 80um at 5um-NMC, other NMC are scaled accordingly
    sys.ly      = 80*ncm.ave(par.iNcm)/ncm.ave(1);        % The box size in y direction is 80um at 5um-NMC, other NMC are scaled accordingly
    sys.lzMn    = 30*ncm.ave(par.iNcm)/ncm.ave(1);        % The minimal length of thickness is  30um at 5um-NMC, other NMC are scaled accordingly
    sys.lzMx    = 120*ncm.ave(par.iNcm)/ncm.ave(1);       % The maximal length of thickness is 120um at 5um-NMC, other NMC are scaled accordingly  
    sys.lz      = sys.lzMn*(sys.lzMx/sys.lzMn)^((ii-1)/(par.nMdl-1));           % Thickness increases as a geometric sequence
    
elseif par.iSty == 3                                      % For massratio study, box size and loading are fixed
    sys.masRt   = (40-ii)/(55+ii);                        % The weight ratio of LPS to NCM
% Scale the system equivalently according to the NMC particle size used in the current study
    sys.load    = 8*10.0*ncm.ave(par.iNcm)/ncm.ave(1);   % The loading is 10mg/cm^2 at 5um-NMC, other NMC are scaled accordingly
    sys.lx      = 200*ncm.ave(par.iNcm)/ncm.ave(1);        % The box size in x direction is 80um at 5um-NMC, other NMC are scaled accordingly
    sys.ly      = 200*ncm.ave(par.iNcm)/ncm.ave(1);        % The box size in y direction is 80um at 5um-NMC, other NMC are scaled accordingly   
    sys.lz      = 100*ncm.ave(par.iNcm)/ncm.ave(1);        % The initial thickness is 80um at 5um-NMC, other NMC are scaled accordingly
    sys.prosty  = 1 - sys.load/(sys.lz*(1+sys.masRt))*(1/ncm.den + sys.masRt/lps.den); % The initial porosity
end

sys.vol     = sys.lx * sys.ly * sys.lz;                                  % The volume of the system
sys.volRt   = sys.masRt * ncm.den/lps.den;                               % The volume ratio of LPS to NCM
lps.vol     = (1 - sys.prosty) * sys.vol * sys.volRt/(1 + sys.volRt);    % The Designed volume of all LPS particles
ncm.vol     = (1 - sys.prosty) * sys.vol * 1/(1 + sys.volRt);            % The Designed volume of all NCM particles

%% 2. Compute actual particle number and volume
% For NMC particles
for itmp = 1 : par.nNcm
%     ncm.s1.nTot = round(6/pi * ncm.vol/ncm.s1.dia^3);
%     ncm.s1.rVol = pi/6 * ncm.s1.dia^3 * ncm.s1.nTot;
    tmp = strcat('ncm.s',num2str(itmp),'.nTot  = round(6/pi * ncm.vol/ncm.s',num2str(itmp),'.dia^3 );');   
    eval(tmp);
    tmp = strcat('ncm.s',num2str(itmp),'.rVol  = pi/6 * ncm.s',num2str(itmp),'.dia^3 * ncm.s',num2str(itmp),'.nTot;');   
    eval(tmp);
end
% For LPS particles
for itmp = 1 : par.nLps
    tmp = strcat('lps.s',num2str(itmp),'.nTot  = round(6/pi * lps.vol/lps.s',num2str(itmp),'.dia^3 );');   
    eval(tmp);
    tmp = strcat('lps.s',num2str(itmp),'.rVol  = pi/6 * lps.s',num2str(itmp),'.dia^3 * lps.s',num2str(itmp),'.nTot;');   
    eval(tmp);
end

