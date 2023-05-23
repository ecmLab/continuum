%% This script is designed for generating models of Ag Nanoparticles under different mass ratios
% "Unit = micro" unit system is used to be consistent as the LAMMPS input file
% Mass: picograms, Distance: micrometers, time:microseconds, density:picograms/micrometer^3

function AgDia = generateAg(cellInf,parInf, iTyp)
%% Model information, including system parameters and material properties of all particles
% cellInf vector: [thickness, xy width]
% parInf  vector: [volRatio, density, particleNumber, aveDia, sgmDia, minDia, maxDia]
    Ag         = struct;
    Ag.volR    = parInf(1);                   % The volume ratio of Ag material
    Ag.den     = parInf(2);                   % The density of Ag material, in unit picogram/micron^3
    Ag.nD      = parInf(3);                   % The number of sampling diameters of Ag particles
    Ag.ave     = parInf(4);                   % The designed average diameter of all Ag, in unit um
    Ag.sgm     = parInf(5);                   % The designed deviation of diameter of all Ag
    Ag.min     = parInf(6);                   % The designed minimal diameter of all Ag, in unit um
    Ag.max     = parInf(7);                   % The designed maximal diameter of all Ag, in unit um
    Ag         = cmp_prob(Ag, Ag.nD, 1);      % Compute the particle distribution

%% Generate modles  
    sys         = struct;
    sys.lz      = cellInf(1);                 % The box size in z direction
    sys.lx      = cellInf(2);                 % The box size in x direction
    sys.ly      = cellInf(2);                 % The box size in y direction
    sys.vol     = sys.lx * sys.ly * sys.lz;   % The volume of the system
    Ag.vol      = Ag.volR * sys.vol;          % The Designed volume of all Ag particles
    Ag          = cmp_nmbr(Ag, Ag.vol);       % Compute actual particle probabilities of all Ag particles

%% Output all particle diameter
    AgDia = zeros(sum(Ag.nDia),1);
    idx    = find(Ag.nDia > 0);
    tmp    = 0;
    for iD = 1 : length(idx)
        AgDia(tmp+1:tmp+Ag.nDia(idx(iD)))     =  Ag.dia(idx(iD))*ones(Ag.nDia(idx(iD)),1);
        tmp = tmp + Ag.nDia(idx(iD));
    end
    
%% Plot particle size distribution when debug mode
    if iTyp == 0
%         figure(1)
%         plot(Ag.dia*1000,100*Ag.pVol,'--k');
%         hold on
        plot(Ag.dia*1000,100*Ag.pVol,'--k', Ag.dia*1000, 100*Ag.dVol,'b');
        legend('Theoretical distribution','Actual distribution')
        xlabel('Diameter (nm)');
        ylabel('Ag size distribution');
    end