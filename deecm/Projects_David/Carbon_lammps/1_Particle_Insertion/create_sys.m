%% This script is designed for generating models for distribution of both LPS and NCM particles
%%% Function: mat_Inpt: define geometry and different types of material and their properties
%%% By: Howard Tu
function [sys,ncmS,ncmL,lpsS,lpsL,carb] = create_sys(ii,c, ncmS, ncmL, lpsS, lpsL, carb)
%% 1. Define System Parameters in sys structure
sys         = struct;                                                                                   % For massratio study, box size and loading are fixed
sys.masRt   = (45-ii)/(55+ii);                                                                          % The weight ratio of LPS and Carbon to NCM
sys.masAt   = 1/(((sys.masRt)/(c*(1+sys.masRt)))-1);                                                    % The weight ratio of Carbon to LPS
sys.load    = 8*10.0;                                                                                   % The loading is 10mg/cm^2
sys.lx      = 50;                                                                                       % The box size in x direction, in unit um 
sys.ly      = 50;                                                                                       % The box size in y direction, in unit um
sys.lz      = 100;                                                                                      % The box size in z direction, in unit um
sys.prosty  = 1 - sys.load/(sys.lz*(1+sys.masRt))*(1/ncmS.den + sys.masRt*((carb.den+lpsS.den*sys.masAt)/((1+sys.masAt)*carb.den*lpsS.den)));               % The initial porosity
sys.vol     = sys.lx * sys.ly * sys.lz;                                                                 % The volume of the system
sys.volRt   = sys.masRt * ncmS.den*(carb.den+sys.masAt*lpsS.den)/(carb.den*lpsS.den*(1+sys.masAt));     % The volume ratio of LPS and Carbon to NCM
sys.volBt   = sys.masAt * lpsS.den/carb.den;                                                            % The volume ratio of Carbon to LPS
ncmS.vol    = (1 - sys.prosty) * sys.vol * 1/(1 + sys.volRt) * ncmS.rt;                                 % The Designed volume of small NCM particles
ncmL.vol    = (1 - sys.prosty) * sys.vol * 1/(1 + sys.volRt) * ncmL.rt;                                 % The Designed volume of large NCM particles
lpsS.vol    = (1 - sys.prosty) * sys.vol * sys.volRt/((1 + sys.volRt)*(1+sys.volBt)) * lpsS.rt;         % The Designed volume of small LPS particles
lpsL.vol    = (1 - sys.prosty) * sys.vol * sys.volRt/((1 + sys.volRt)*(1+sys.volBt)) * lpsL.rt;         % The Designed volume of large LPS particles
carb.vol    = (1 - sys.prosty) * sys.vol * sys.volRt*sys.volBt/((1 + sys.volRt)*(1+sys.volBt));         % The Designed volume of carbon particles

%% 2. Compute actual particle number and volume
% For small NMC particles
ncmS.nTot = round(6/pi * ncmS.vol/ncmS.dia^3);  % The number of NCM small particles
ncmS.rVol = pi/6 * ncmS.dia^3 * ncmS.nTot;      % The real volume of all small NCM particles

% For large NMC particles
ncmL.nTot = round(6/pi * ncmL.vol/ncmL.dia^3);  % The number of NCM large particles
ncmL.rVol = pi/6 * ncmL.dia^3 * ncmL.nTot;      % The real volume of all large NCM particles

% For small LPS particles
lpsS.nTot = round(6/pi * lpsS.vol/lpsS.dia^3);  % The number of LPS small particles
lpsS.rVol = pi/6 * lpsS.dia^3 * lpsS.nTot;      % The real volume of all small LPS particles

% For large LPS particles
lpsL.nTot = round(6/pi * lpsL.vol/lpsL.dia^3);  % The number of LPS large particles
lpsL.rVol = pi/6 * lpsL.dia^3 * lpsL.nTot;      % The real volume of all large LPS particles

% For carbon particles
carb.nTot = round(6/pi * carb.vol/carb.dia^3);  % The number of carbon particles
carb.rVol = pi/6 * carb.dia^3 * carb.nTot;      % The real volume of all carbon particles

% Other parameters
sys.mass = (lpsS.rVol + lpsL.rVol)*lpsS.den + (ncmS.rVol+ncmL.rVol)*ncmS.den + carb.rVol*carb.den;     % The total mass of the model, in unit pico-gram
