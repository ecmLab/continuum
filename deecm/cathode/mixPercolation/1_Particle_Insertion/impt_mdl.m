%% This script is designed for generating models for distribution of both LPS and NCM particles
%%% Function: impt_mat: define geometry and material properties of all particles
%%% By: Howard Tu

%% 1.Define properties of small NCM particles in ncm structure
ncmS         = struct;
ncmS.nD      = 20;            % The number of sampling diameters of NCM particles
ncmS.den     = 4.85;          % The density of NCM material, in unit picogram/micron^3
ncmS.typ     = 1;             % Label NCM particles as Type 1 particles

%% 2.Define properties of large NCM particles in ncm structure
ncmL         = struct;
ncmL.nD      = 20;            % The number of sampling diameters of NCM particles
ncmL.den     = 4.85;          % The density of NCM material, in unit picogram/micron^3
ncmL.typ     = 1;             % Label NCM particles as Type 1 particles

%% 3.Define properties of small LPS particles in lpsS structure
lpsS        = struct;
lpsS.nD     = 20;            % The number of sampling diameters of small LPS particles
lpsS.den    = 1.87;          % The density of LPS material, in unit picogram/micron^3
lpsS.typ    = 2;             % Label LPS particles as Type 2 particles

%% 4.Define properties of large LPS particles in lpsS structure
lpsL        = struct;
lpsL.nD     = 20;            % The number of sampling diameters of small LPS particles
lpsL.den    = 1.87;          % The density of LPS material, in unit picogram/micron^3
lpsL.typ    = 2;             % Label LPS particles as Type 2 particles

%% 5.Define properties of carbon particles in carbon structure
carb        = struct;
carb.nD      = 20;            % The number of sampling diameters of large LPS particles
carb.den     = 2.26;          % The density of LPS material, in unit picogram/micron^3
carb.typ     = 3;             % Label carbon particles as Type 3 particles
carb.dia     = 0.25;           % Diameter of carbon particles

%% 6.Define ghost particles at the top of the system to apply pressure
gst         = struct;
gst.den     = 10;            % The density of ghost material
gst.typ     = 4;

