%% This script is designed for generating models of both LPS and NCM particles under different mass ratios of LPS to NCM
%%% By: Howard Tu
clc; clear;
% "Unit = micro" unit system is used to be consistent as the LAMMPS input file
%  Mass: picograms, Distance: micrometers, time:microseconds, density:picograms/micrometer^3

%% Define variables
nWt    = 1;                    % The number of different mass ratios (LPS:NCM) created, following the equation: (40-ii)/(55+ii), ii=1 to nWt

%% Import model information, including system parameters and material properties of all particles
impt_mdl;
for nc=1:1
sm=[1 1 2 2];
lm=[3 4 3 4];
lpsS.dia     = sm(nc);            % The designed diameter of small LPS particles, in unit um
lpsL.dia     = lm(nc);            % The designed diameter of large LPS particles, in unit um

for ratio=1:1
    A=[0.2 0.4 0.5 0.6 0.8];
    lpsS.rt=A(ratio);
    lpsL.rt=1-A(ratio);
    ncmS.rt=0.5;
    ncmL.rt=0.5;
    lambda=2;
    meanlps=(lpsS.rt*lpsS.dia+lpsL.rt*lpsL.dia);
    ncmS.dia=lambda*meanlps;
    ncmL.dia=lambda*meanlps;
%% Generate nMdl data files for different mass ratios
cc=1;
c=0.01*cc; %Percentage of carbon particles
for iwt = 1:nWt    
% 0. Generate folder of different files
    fid = strcat('../2LPS_',num2str(lpsS.dia),'_',num2str(lpsL.dia),'/Ratio',num2str(ratio),'/massratio/mr',num2str(iwt),'/');
    stmp = ['mkdir ' fid];
    eval(stmp);
    
% 1. Generate system, including box and system level variables
    [sys,ncmS,ncmL,lpsS,lpsL,carb] = create_sys(iwt,c, ncmS, ncmL, lpsS, lpsL, carb);
    
% 2. Insert all Particles into the box, stored in sys.cord matrix:
    % 1st colume: id of particle;
    % 2th column: particle type: LPS or NCM;
    % 3th column: diameter of particle;
    % 4th column: density of particle
    % 5-7th, coordinate of each particle
    % Note: it would be more efficient if bigger particles are inserted first, then smaller particles
    % Initialize the cordinate matrix with one biggest NCM particle, Particle's boundary should not outside the wall
    if lpsL.dia < ncmL.dia
        sys.cord  = [ 1, ncmL.typ, ncmL.dia, ncmL.den, (sys.lx - ncmL.dia)*rand + ncmL.dia/2,...
            (sys.ly - ncmL.dia)*rand + ncmL.dia/2, (sys.lz - ncmL.dia)*rand + ncmL.dia/2];
        sys.cord  = insrtPtc(sys, ncmL, ncmL.typ, ncmL.den, 1);  % Insert all large NCM particles
        sys.cord  = insrtPtc(sys, ncmS, ncmS.typ, ncmS.den, 1);  % Insert all small NCM particles
        sys.cord  = insrtPtc(sys, lpsL, lpsL.typ, lpsL.den, 1);  % Insert all large LPS particles
        sys.cord  = insrtPtc(sys, lpsS, lpsS.typ, lpsS.den, 1);  % Insert all small LPS particles
        sys.cord  = insrtPtc(sys, carb, carb.typ, carb.den, 1);  % Insert all carbon particles
        gst.max   = carb.dia;                                    % The maximal possible size for ghost particless
    else
        sys.cord  = [ 1, lpsL.typ, lpsL.dia, lpsL.den, (sys.lx - lpsL.dia)*rand + lpsL.dia/2,...
            (sys.ly - lpsL.dia)*rand + lpsL.dia/2, (sys.lz - lpsL.dia)*rand + lpsL.dia/2];
        sys.cord  = insrtPtc(sys, lpsL, lpsL.typ, lpsL.den, 1);  % Insert all large LPS particles
        sys.cord  = insrtPtc(sys, ncmL, ncmL.typ, ncmL.den, 1);  % Insert all large NCM particles
        sys.cord  = insrtPtc(sys, ncmS, ncmS.typ, ncmS.den, 1);  % Insert all small NCM particles
        sys.cord  = insrtPtc(sys, lpsS, lpsS.typ, lpsS.den, 1);  % Insert all small LPS particles
        sys.cord  = insrtPtc(sys, carb, carb.typ, carb.den, 1);  % Insert all carbon particles
        sys.cord(:,2:end)  = sortrows(sys.cord(:,2:end),1);      % Make sure nmc particles always starts first
        gst.max   = carb.dia;                                    % The maximal possible size for ghost particles
    end

% 3. Generate the coordinates of all ghost particles
    gst.nxy  = ceil(3*(sqrt(2)-1) * sys.lx/gst.max);  % The number of ghost particles in X and Y direction for current particle size
    gst.dia  = sys.lx/gst.nxy;                        % The diameter of ghost particles for current particle size
    gst.cord = zeros(gst.nxy^2,7);                    % initialize the ghost particles
    gst.cord(:,1) = [size(sys.cord,1)+1 : size(sys.cord,1)+gst.nxy^2];   % Ghost particles index after all LPS and NCM particles
    gst.cord(:,2) = gst.typ;
    gst.cord(:,3) = gst.dia;
    gst.cord(:,4) = gst.den;
    gst.cord(:,7) = sys.lz + gst.dia/2;               % All ghost particles located at the top of the system with a radius distance
    for ig = 1 : gst.nxy
        for jg = 1 : gst.nxy
            gst.cord((ig-1)*gst.nxy+jg,5) = gst.dia/2 + gst.dia*(jg-1);
            gst.cord((ig-1)*gst.nxy+jg,6) = gst.dia/2 + gst.dia*(ig-1);
        end
    end

% 4. Output all particles in LAMMPS data file format
    stmp = strcat(fid,'mdl.data');
    fileID = fopen(stmp,'w');
    fprintf(fileID,'LAMMPS data file\n');
    fprintf(fileID,'\n');
    fprintf(fileID,'%6d %6s\n',size(sys.cord,1)+size(gst.cord,1),'atoms');
    fprintf(fileID,'\n');
    fprintf(fileID,'   4 atom types\n');
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
    
end   % End of iwt

end
end
