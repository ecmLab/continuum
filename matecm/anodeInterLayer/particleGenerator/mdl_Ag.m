%% This script is designed for generating models of Ag Nanoparticles under different mass ratios
% "Unit = micro" unit system is used to be consistent as the LAMMPS input file
% Mass: picograms, Distance: micrometers, time:microseconds, density:picograms/micrometer^3

%% Define variables
iTyp       = 0;                     % iTyp==0 is for debugging, no real model will be generated; iTyp==1 is for generate model

%% Model information, including system parameters and material properties of all particles
Ag         = struct;
Ag.nD      = 50;            % The number of sampling diameters of Ag particles
Ag.den     = 4.85;          % The density of Ag material, in unit picogram/micron^3
Ag.ave     = 0.04;          % The designed average diameter of all Ag, in unit um
Ag.sgm     = 0.01;          % The designed deviation of diameter of all Ag
Ag.min     = 0.02;          % The designed minimal diameter of all Ag, in unit um
Ag.max     = 0.1;           % The designed maximal diameter of all Ag, in unit um
Ag         = cmp_prob(Ag,Ag.nD,1); % Compute the particle distribution

%% Generate modles  
sys         = struct;
sys.load    = 0.1;       % The loading is 1mg/cm^2 at 5um-NMC, other NMC are scaled accordingly
sys.lx      = 1;                                % The box size in x direction
sys.ly      = 1;                                % The box size in y direction
sys.lz      = 1;                                % The box size in z direction
sys.prosty  = 1 - sys.load/(sys.lz*Ag.den)*10;  % The initial porosity
sys.vol     = sys.lx * sys.ly * sys.lz;         % The volume of the system
Ag.vol      = (1 - sys.prosty) * sys.vol;       % The Designed volume of all Ag particles
Ag          = cmp_nmbr(Ag,Ag.vol); % Compute actual particle probabilities of all Ag particles

%Plot particle size distribution when debug mode
if iTyp == 0
    figure(1)
    plot(Ag.dia,100*Ag.pVol,'--k', Ag.dia, 100*Ag.dVol,'b');
    legend('Theoretical distribution','Actual distribution')
    xlabel('Diameter (\mum)');
    ylabel('Probability');
end