%% This script is designed for generating models of the DRX cathode+Carbon used for Haegyum's paper
%  Different models under different particle distributions of DRX to C
%%% By: Howard Tu, 01/30/2026
clc; clear;
% "Unit = micro" unit system is used to be consistent with the LAMMPS input file
% Mass: picograms, Distance: micrometers, time:microseconds, density:picograms/micrometer^3
nWt    = 1;  

for isz  = 1:10
for iwt = 1:11
for lo = 1:8

    ll=[2 4 8 10 12 14 16 18];
%% Define variables
% Four different Ag size distributions are considered:
% isz = 1 means small size + fine mix: Ag_ave = 40nm, narrow lognorm distribution
% isz = 2 means small size + hand mix: Ag_ave = 40nm, wide lognorm distribution
% isz = 3 means large size + fine mix: Ag_ave = 80nm, narrow lognorm distribution
% isz = 4 means large size + hand mix: Ag_ave = 80nm, wide lognorm distribution
% 1. Generate involved particles with corresponding material properties and geometries
[ncm,lps,gst] = particle_info(isz);
% 2. Generate system, including box and system level variables
[sys,ncm,lps] = create_sys(ncm, lps,iwt,ll(lo));

% 0. Generate folder of different files
    fid = strcat('./Load_',num2str(ll(lo)),'/LPS25_',num2str(lps.sgm),'_NMC5_',num2str(ncm.sgm),'/massratio/mr',num2str(iwt),'/');
    stmp = ['mkdir ' fid];
    eval(stmp);

%fprintf('The total particles are Ag = %6.4f , C = %6.4f\n',ag.nTot, ac.nTot);
%fprintf('The initial porosity is : %6.4f\n',1-(ag.rVol+ac.rVol)/sys.vol);

% 3. Insert all Particles into the box, stored in sys.cord matrix:
  % 1st colume: id of particle;
  % 2th column: particle type: Ag or aC;
  % 3th column: diameter of particle;
  % 4th column: density of particle
  % 5-7th, coordinate of each particle
  % Note: it would be more efficient if bigger particles are inserted first, then smaller particles
  % Initialize the cordinate matrix with one biggest Ag particle, Particle's boundary should not outside the wall
if lps.ave < ncm.ave
    sys.cord  = [1, ncm.typ, ncm.dia(end), ncm.den, (sys.lx - ncm.dia(end))*rand + ncm.dia(end)/2,...
                 (sys.ly - ncm.dia(end))*rand + ncm.dia(end)/2, (sys.lz - ncm.dia(end))*rand + ncm.dia(end)/2];
    sys.cord  = insrtPtc(sys, ncm, ncm.typ, ncm.den, 1);     % Insert all Ag particles
    sys.cord  = insrtPtc(sys, lps, lps.typ, lps.den, 1);     % Insert all C particles
    gst.max   = min(lps.dia);                                     % The maximal possible size for ghost particles
else
    sys.cord  = [1, lps.typ, lps.dia(end), lps.den, (sys.lx - lps.dia(end))*rand + lps.dia(end)/2,...
                 (sys.ly - lps.dia(end))*rand + lps.dia(end)/2, (sys.lz - lps.dia(end))*rand + lps.dia(end)/2];
    sys.cord  = insrtPtc(sys, lps, lps.typ, lps.den, 1);  % Insert all C particles
    sys.cord  = insrtPtc(sys, ncm, ncm.typ, ncm.den, 1);  % Insert all Ag particles
    sys.cord(:,2:end)  = sortrows(sys.cord(:,2:end),1);      % Make sure Ag particles are always listed before C particles in the input data first
    gst.max   = min(ncm.dia);                                  % The maximal possible size for ghost particles
end

% 4. Generate the coordinates of all ghost particles
gst.nxy  = ceil(3*(sqrt(2)-1) * sys.lx/gst.max);  % The number of ghost particles in X and Y direction for current particle size
gst.dia  = sys.lx/gst.nxy;                      % The diameter of ghost particles for current particle size
gst.cord = zeros(gst.nxy^2,7);                  % initialize the ghost particles
gst.cord(:,1) = [size(sys.cord,1)+1 : size(sys.cord,1)+gst.nxy^2];   % Ghost particles index after all Ag and C particles
gst.cord(:,2) = gst.typ;
gst.cord(:,3) = gst.dia;
gst.cord(:,4) = gst.den;
gst.cord(:,7) = sys.lz + gst.dia/2;              % All ghost particles located at the top of the system with a radius distance
for ig = 1 : gst.nxy
    for jg = 1 : gst.nxy
        gst.cord((ig-1)*gst.nxy+jg,5) = gst.dia/2 + gst.dia*(jg-1);
        gst.cord((ig-1)*gst.nxy+jg,6) = gst.dia/2 + gst.dia*(ig-1);
    end
end

% 5. Output all particles in LAMMPS data file format
stmp = strcat(fid,'mdl.data');
fileID = fopen(stmp,'w');
fprintf(fileID,'LAMMPS data file\n');
fprintf(fileID,'\n');
fprintf(fileID,'%6d %6s\n',size(sys.cord,1)+size(gst.cord,1),'atoms');
fprintf(fileID,'\n');
fprintf(fileID,'   3 atom types\n');
fprintf(fileID,'\n');
fprintf(fileID,'%23.19f %23.19f %3s %3s\n',0.0, sys.lx, 'xlo', 'xhi');
fprintf(fileID,'%23.19f %23.19f %3s %3s\n',0.0, sys.ly, 'ylo', 'yhi');
fprintf(fileID,'%23.19f %23.19f %3s %3s\n',0.0, sys.lz+gst.dia, 'zlo', 'zhi');
fprintf(fileID,'\n');
fprintf(fileID,'Atoms # sphere\n');
fprintf(fileID,'\n');
fprintf(fileID,'%5d %4d %8.5f %8.5f %20.15f %20.15f %20.15f\n',sys.cord');
fprintf(fileID,'%5d %4d %8.5f %8.5f %20.15f %20.15f %20.15f\n',gst.cord');
fclose(fileID);
end
end
end