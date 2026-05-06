%% This script is designed for generating microstructure of SE with porosity with different density and pore distribution
%%% This is the first version created on 11/15/2019
%%% By: Howard Tu
clc; clear;
% "Unit = micro" unit system is used to be consistent as the LAMMPS input file
% Mass: picograms, Distance: micrometers, time:microseconds, density:picograms/micrometer^3

%% Define variables
% Define parameters in par structure
par         = struct;
par.iSty    = 3;                     % Study types: 1 means box size study, 2 means loading study, 3 means massratio study
par.iTyp    = 1;                     % iTyp==0 is for debugging, no real model will be generated; iTyp==1 is for generate model
par.icnf    = 1;                     % current configurations

% Define System Parameters in sys structure
sys         = struct;
sys.prosty  = 0.02;                                    % The porosity of the SE
sys.lx      = 100;                                    % The box size in x direction
sys.ly      = 40;                                     % The box size in y direction
sys.lz      = 40;                                     % The box size in z direction
sys.vol     = sys.lx * sys.ly * sys.ly;                % The volume of the system
sys.lm      = 1;                                       % The Characteristic length for the mesh

% Define Pore parameters in por structure
por         = struct;
por.vol     = sys.prosty * sys.vol;                   % The Designed volume of all pore particles
por.min     = 0.5;                                    % The minimum radius of pore particles, in unit micron
por.max     = 0.5;                                   % The maximum radius of pore particles, in unit micron
por.nD      = 10;                                     % The number of sampling radius of pore particles
tmp         = linspace(por.min,por.max,por.nD+1);
por.rds     = tmp(1:end-1) + (por.max-por.min)/(2*por.nD);  % The sampling diameters of all pores
por.nRds    = floor(3/(4*pi)* por.vol / sum(por.rds.^3))*ones(1,por.nD);       % Approximate number of pore particles for each diameter


% 2. Generate models
if par.iTyp == 1       % Only generate model when iTyp==1 to save time
    %% Insert all pore Particles into the box, stored in sys.cord matrix:
    % 1st colume: id of particle;
    % 2th column: diameter of particle;
    % 3-4th, coordinate of each particle
    % Note: it would be more efficient if bigger particles are inserted first, then smaller particles
    % Initialize the cordinate matrix with one biggest particle, Particle's boundary should not outside the wall
    sys.cord  = [ 1, por.rds(end), (sys.lx - 4*por.rds(end))*rand + 2*por.rds(end), (sys.ly - 4*por.rds(end))*rand + 2*por.rds(end), (sys.lz - 4*por.rds(end))*rand + 2*por.rds(end)];
    sys.cord  = insrtPtc(sys, por);     % Insert all other pore particles
    
    %% Output all particles in GMesh data file format
    % Output filename
    fout = fopen('../mdl1.geo','w');
    
    % Write parameters
    str = strcat('m = ',num2str(sys.lm),'; // mesh characteristic length');
    fprintf(fout,'%s\n',str);
    fprintf(fout,'%s\n','SetFactory("OpenCASCADE");');
    
    % Write the SE box
    fprintf(fout,'\n');
    str = strcat('Box(1) = {0,0,0,',num2str(sys.lx),',',num2str(sys.ly),',',num2str(sys.lz),'};');
    fprintf(fout,'%s\n',str);
    
    % Write the coordinate of pore center and radius
    fprintf(fout,'\n');
    for i = 1 : length(sys.cord)
        % Create the coordinate of pore center and boundary points
        str = strcat('Sphere(',num2str(i+1),') = {',num2str(sys.cord(i,3), '%10.6f'),',',num2str(sys.cord(i,4), '%10.6f'),',',num2str(sys.cord(i,5), '%10.6f'),',',num2str(sys.cord(i,2), '%10.6f'),'};');
        fprintf(fout,'%s\n',str);
      
    end
    
%     % Boolean Union to all pores
%     fprintf(fout,'\n');
%     str = strcat('BooleanUnion(',num2str(i+2),') = {Volume{2:',num2str(i+1),'};Delete;}{Volume{2:',num2str(i+1),'};Delete;};');
%     fprintf(fout,'%s\n',str);

    % Boolean Difference: Substract all pores from the SE
    fprintf(fout,'\n');
    str = strcat('BooleanDifference(',num2str(i+2),') = {Volume{1};Delete;}{Volume{2:',num2str(i+1),'};Delete;};');
    fprintf(fout,'%s\n',str);    
    
    % Define characteristic length
    fprintf(fout,'\n');
    str = strcat('Characteristic Length{ PointsOf{ Volume{',num2str(i+2),'}; } } = ',num2str(sys.lm),';');
    fprintf(fout,'%s\n',str);
      
    % Create the Physical names
    fprintf(fout,'\n');
    str = strcat('Physical Surface("left") = {1};');
    fprintf(fout,'%s\n',str);
    
    str = strcat('Physical Surface("right") = {6};');
    fprintf(fout,'%s\n',str);
    
    fprintf(fout,'%s\n','nS = newreg;');
    str = strcat('Physical Surface("pore") = {7:nS-1};');
    fprintf(fout,'%s\n',str);

    str = strcat('Physical Volume("blockSE") = {',num2str(i+2),'};');
    fprintf(fout,'%s\n',str);
    
    fclose(fout);
    
end    % End of iType

