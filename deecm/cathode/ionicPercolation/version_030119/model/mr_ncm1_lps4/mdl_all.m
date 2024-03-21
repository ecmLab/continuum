%% This script is designed for generating models of both LPS and NCM particles under different mass ratios
%%% This is the third version created on April,16st,2018, primary for the Group Talk on April, 20
%%% By: Howard Tu
clc; clear;
% "Unit = micro" unit system is used to be consistent as the LAMMPS input file
% Mass: picograms, Distance: micrometers, time:microseconds, density:picograms/micrometer^3

%% Define variables
% Define parameters in par structure
par         = struct;
par.nNcm    = 3;                     % Number of totoal NMC sizes
par.nLps    = 9;                     % Number of totoal LPS sizes
par.iSty    = 3;                     % Study types: 1 means box size study, 2 means loading study, 3 means massratio study
par.iTyp    = 0;                     % iTyp==0 is for debugging, no real model will be generated; iTyp==1 is for generate model
par.iNcm    = 1;                     % The size of NCM used, 1 means 5um, 2 means 10um, 3 means 12um
par.icnf    = 1;                     % current configurations
if par.iSty == 1
   par.nMdl  = 11;                   % 11 models are generated for boxsize study
   par.ncmMr = 80;                   % Mass ratio of NCM of current study
elseif par.iSty == 2
   par.nMdl  = 20;                   % 20 models are generated for loading study
   par.ncmMr = 80;                   % Mass ratio of NCM of current study
   par.prsty = 0.7;                  % Designed initial porosity for current study   
elseif par.iSty == 3 
   par.nMdl = 35;                    % 35 models are generated for massratio study
end

% Define results in rst structure
rst     = struct;
rst.ldg = zeros(par.nMdl,par.nLps);         % Cathode loading, in unit mg/cm^2
rst.prs = zeros(par.nMdl,par.nLps);         % Porosity
rst.mrt = zeros(par.nMdl,par.nLps+1);       % Mass ratio of LPS:NCM
rst.nrt = zeros(par.nMdl,par.nLps);         % Number ratio of LPS:NCM
rst.num = zeros(par.nMdl,par.nLps);         % Total number of particles
if par.iSty == 1
   rst.mrt(:,par.nLps+1) = (95-par.ncmMr)/par.ncmMr;   % Sotre designed mass ratio of LPS:NCM
elseif par.iSty == 2
   rst.mrt(:,par.nLps+1) = (95-par.ncmMr)/par.ncmMr;   % Sotre designed mass ratio of LPS:NCM
elseif par.iSty == 3
   rst.mrt(:,par.nLps+1) = (55+[1:par.nMdl])./(40-[1:par.nMdl]);   % Sotre designed mass ratio of LPS:NCM
end

%% Import model information, including system parameters and material properties of all particles
[ncm,lps,gst] = impt_mdl(par);

%% Generate nMdl modles
for imdl = 1:par.nMdl

% 0. Generate folder of different models
    if     par.iSty == 1                  % For Boxsize convergence study
        fid = strcat('../../',num2str(par.icnf),'conf/data/box/bx',num2str(imdl),'/');
    elseif par.iSty == 2                  % For loading study
        fid = strcat('../../',num2str(par.icnf),'conf/data/loading/ld',num2str(imdl),'/');
    elseif par.iSty == 3                  % For massratio study
        fid = strcat('../../',num2str(par.icnf),'conf/data/massRatio/mr',num2str(imdl),'/');
    end

    if par.iTyp == 1
        stmp = ['mkdir ' fid];
        eval(stmp);
    end   

% 1. Generate system, including box and system level variables
    [sys,ncm,lps] = create_sys(imdl, par, ncm, lps);
   % The NMC particle used in the current study
    stmp          = strcat('ncm.sz = ncm.s',num2str(par.iNcm),';');
    eval(stmp); 

% 2. Generate models of different LPS sizes of current mass ratio and current NMC particles
    for isz = 1 : 1
        stmp     = strcat('lps.sz = lps.s',num2str(isz),';');
        eval(stmp);                                                        % Dealing with particle size isz
        sys.mass = lps.sz.rVol*lps.den + ncm.sz.rVol*ncm.den;              % The total mass of the model, in unit pico-gram
        rst.mrt(imdl,isz) = (lps.sz.rVol*lps.den)/(ncm.sz.rVol*ncm.den);   % Sotre real mass ratio of LPS:NCM
        rst.prs(imdl,isz) = 1 - (lps.sz.rVol+ncm.sz.rVol)/sys.vol;         % Sotre real porosity
        rst.ldg(imdl,isz) = sys.mass/(sys.lx*sys.ly)*0.1;                  % Sotre real cathod loading
        rst.num(imdl,isz) = sum(ncm.sz.nDia);           % Sotre total number of particles
        rst.nrt(imdl,isz) = sum(ncm.sz.nDia)/sum(lps.sz.nDia);             % Sotre total number of particles
                
%         fprintf('The real Mass ratio of LPS to NCM is : %6.4f\n',rst.mrt(itk,isz));
%         fprintf('The real Porosity is : %6.4f\n',rst.prs(itk,isz));
%         fprintf('The real Mass per area is : %6.4f\n',rst.ldg(imdl,isz));
        fprintf('\n\n Dealing with model : %6d of LPS size: %6d\n',imdl, isz);

       if par.iTyp == 1       % Only generate model when iTyp==1 to save time
    %% 6. Insert all Particles into the box, stored in sys.cord matrix:
        % 1st colume: id of particle;
        % 2th column: particle type: LPS or NCM;
        % 3th column: diameter of particle;
        % 4th column: density of particle
        % 5-7th, coordinate of each particle
        % Note: it would be more efficient if bigger particles are inserted first, then smaller particles
        % Initialize the cordinate matrix with one biggest NCM particle, Particle's boundary should not outside the wall
        if lps.sz.ave < ncm.sz.ave
            sys.cord  = [ 1, ncm.typ, ncm.sz.dia(end), ncm.den, (sys.lx - ncm.sz.dia(end))*rand + ncm.sz.dia(end)/2,...
                          (sys.ly - ncm.sz.dia(end))*rand + ncm.sz.dia(end)/2, (sys.lz - ncm.sz.dia(end))*rand + ncm.sz.dia(end)/2];
            sys.cord  = insrtPtc(sys, ncm.sz, ncm.typ, ncm.den, 1);     % Insert all NCM particles
            sys.cord  = insrtPtc(sys, lps.sz, lps.typ, lps.den, 1);     % Insert all LPS particles
            gst.max   = lps.sz.min;                                     % The maximal possible size for ghost particles
        else
            sys.cord  = [ 1, lps.typ, lps.sz.dia(end), lps.den, (sys.lx - lps.sz.dia(end))*rand + lps.sz.dia(end)/2,...
                          (sys.ly - lps.sz.dia(end))*rand + lps.sz.dia(end)/2, (sys.lz - lps.sz.dia(end))*rand + lps.sz.dia(end)/2];
            sys.cord  = insrtPtc(sys, lps.sz, lps.typ, lps.den, 1);  % Insert all LPS particles
            sys.cord  = insrtPtc(sys, ncm.sz, ncm.typ, ncm.den, 1);  % Insert all NCM particles
            sys.cord(:,2:end)  = sortrows(sys.cord(:,2:end),1);      % Make sure nmc particles always starts first
            gst.max   = ncm.sz.min;                                  % The maximal possible size for ghost particles
        end

    %% 7. Generate the coordinates of all ghost particles
        gst.nxy  = ceil(3*(sqrt(2)-1) * sys.lx/gst.max);  % The number of ghost particles in X and Y direction for current particle size
        gst.dia  = sys.lx/gst.nxy;                      % The diameter of ghost particles for current particle size
        gst.cord = zeros(gst.nxy^2,7);                  % initialize the ghost particles
        gst.cord(:,1) = [size(sys.cord,1)+1 : size(sys.cord,1)+gst.nxy^2];   % Ghost particles index after all LPS and NCM particles
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
        
    %% 8. Output all particles in LAMMPS data file format
        if par.iSty < 3            % For boxsize and loading study
            stmp = strcat(fid,'mdl_ncm',num2str(par.iNcm),'_lps',num2str(isz),'_mr',num2str(par.ncmMr),'.data');
        elseif par.iSty == 3       % For massratio study 
            stmp = strcat(fid,'mdl_ncm',num2str(par.iNcm),'_lps',num2str(isz),'.data');
        end
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

       end    % End of iType
         
    end    % End of isz
    
end   % End of imdl

% fileID = fopen('../5conf/data/thickness/massRatio.txt','w');
% fprintf(fileID,'%10.6f %10.6f %10.6f %10.6f %10.6f\n',rst.mrt');
% fclose(fileID);
