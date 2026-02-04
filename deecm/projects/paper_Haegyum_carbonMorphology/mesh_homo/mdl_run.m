%% This script is designed for generating models of the DRX cathode + Carbon used for Haegyum's project
%  Different models under different particle distributions of DRX to carbon
%%% By: Howard Tu, 01/30/2026
clc; clear;
% "Unit = micro" unit system is used to be consistent with the LAMMPS input file
% Mass: picograms, Distance: micrometers, time:microseconds, density:picograms/micrometer^3
nWt    = 11;  

for iwt = 1:nWt

%% Define variables
[drx,carbon,gst] = particle_info();
% 2. Generate system, including box and system level variables
[sys,drx,carbon] = create_sys(drx, carbon,iwt);

% 0. Generate folder of different files
    fid = strcat('massratio/mr',num2str(iwt),'/');
    stmp = ['mkdir ' fid];
    eval(stmp);

fprintf('The total particles are DRX = %6.4f , C = %6.4f\n',drx.nTot, carbon.nTot);
fprintf('The initial porosity is : %6.4f\n',1-(drx.rVol+carbon.rVol)/sys.vol);

% 3. Insert all Particles into the box, stored in sys.cord matrix:
  % 1st colume: id of particle;
  % 2th column: particle type: Ag or aC;
  % 3th column: diameter of particle;
  % 4th column: density of particle
  % 5-7th, coordinate of each particle
  % Note: it would be more efficient if bigger particles are inserted first, then smaller particles
  % Initialize the cordinate matrix with one biggest Ag particle, Particle's boundary should not outside the wall
if carbon.ave < drx.ave
    sys.cord  = [1, drx.typ, drx.dia(end), drx.den, (sys.lx - drx.dia(end))*rand + drx.dia(end)/2,...
                 (sys.ly - drx.dia(end))*rand + drx.dia(end)/2, (sys.lz - drx.dia(end))*rand + drx.dia(end)/2];
    sys.cord  = insrtPtc(sys, drx, drx.typ, drx.den, 1);     % Insert all DRX particles
    sys.cord  = insrtPtc(sys, carbon, carbon.typ, carbon.den, 1);     % Insert all carbon particles
    gst.max   = min(carbon.dia);                                     % The maximal possible size for ghost particles
else
    sys.cord  = [1, carbon.typ, carbon.dia(end), carbon.den, (sys.lx - carbon.dia(end))*rand + carbon.dia(end)/2,...
                 (sys.ly - carbon.dia(end))*rand + carbon.dia(end)/2, (sys.lz - carbon.dia(end))*rand + carbon.dia(end)/2];
    sys.cord  = insrtPtc(sys, carbon, carbon.typ, carbon.den, 1);  % Insert all carbon particles
    sys.cord  = insrtPtc(sys, drx, drx.typ, drx.den, 1);           % Insert all DRX particles
    sys.cord(:,2:end)  = sortrows(sys.cord(:,2:end),1);            % Make sure DRX particles are listed before carbon in the input data
    gst.max   = min(drx.dia);                                      % The maximal possible size for ghost particles
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