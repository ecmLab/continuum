%% This script is designed for generating models of the AgC interlayer used for YaoYan's paper
%  Different models under different particle distributions of Ag to C
%%% By: Howard Tu, 06/10/2023
clc; clear;
% "Unit = micro" unit system is used to be consistent with the LAMMPS input file
% Mass: picograms, Distance: micrometers, time:microseconds, density:picograms/micrometer^3

%% Define variables
% Four different Ag size distributions are considered:
% isz = 1 means small size + fine mix: Ag_ave = 40nm, narrow lognorm distribution
% isz = 2 means small size + hand mix: Ag_ave = 40nm, wide lognorm distribution
% isz = 3 means large size + fine mix: Ag_ave = 80nm, narrow lognorm distribution
% isz = 4 means large size + hand mix: Ag_ave = 80nm, wide lognorm distribution
% 1. Generate involved particles with corresponding material properties and geometries

dev=[0.1 0.3 0.5 0.7 0.9 1 1.5 2];

for i=1:length(dev)

[lps,gst] = particle_info(dev(i));
% 2. Generate system, including box and system level variables
[sys,lps] = create_sys(lps);



% 0. Generate folder of different files
    fid = strcat('../LPS/D3_dev',num2str(dev(i)),'/');
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
x1=(sys.lx/2 - lps.dia(end))*randi([-1 1],1)*rand + lps.dia(end)/2;
y1=(sqrt((sys.ly/2)^2 - x1^2)-lps.dia(end))*randi([-1 1],1)*rand + lps.dia(end)/2;
sys.cord  = [1, lps.typ, lps.dia(end), lps.den, x1, y1, (sys.lz - lps.dia(end))*rand + lps.dia(end)/2];
sys.cord  = insrtPtc(sys, lps, lps.typ, lps.den, 1);  % Insert all C particles
sys.cord(:,2:end)  = sortrows(sys.cord(:,2:end),1);      % Make sure Ag particles are always listed before C particles in the input data first
gst.max   = min(lps.dia);                                  % The maximal possible size for ghost particles


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
        gst.cord((ig-1)*gst.nxy+jg,5) = -sys.lx/2+gst.dia/2 + gst.dia*(jg-1);
        gst.cord((ig-1)*gst.nxy+jg,6) = -sys.ly/2+gst.dia/2 + gst.dia*(ig-1);
    end
end

% 5. Output all particles in LAMMPS data file format
stmp = strcat(fid,'mdl.data');
fileID = fopen(stmp,'w');
fprintf(fileID,'LAMMPS data file\n');
fprintf(fileID,'\n');
fprintf(fileID,'%6d %6s\n',size(sys.cord,1)+size(gst.cord,1),'atoms');
fprintf(fileID,'\n');
fprintf(fileID,'   2 atom types\n');
fprintf(fileID,'\n');
fprintf(fileID,'%23.19f %23.19f %3s %3s\n',-sys.lx/2, sys.lx/2, 'xlo', 'xhi');
fprintf(fileID,'%23.19f %23.19f %3s %3s\n',-sys.ly/2, sys.ly/2, 'ylo', 'yhi');
fprintf(fileID,'%23.19f %23.19f %3s %3s\n',0.0, sys.lz+gst.dia, 'zlo', 'zhi');
fprintf(fileID,'\n');
fprintf(fileID,'Atoms # sphere\n');
fprintf(fileID,'\n');
fprintf(fileID,'%5d %4d %8.5f %8.5f %20.15f %20.15f %20.15f\n',sys.cord');
fprintf(fileID,'%5d %4d %8.5f %8.5f %20.15f %20.15f %20.15f\n',gst.cord');
fclose(fileID);
end