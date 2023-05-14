%% This code is developed for the Li-dendrite experimental paper
%   Howard Tu, 02/19/2022
clc; clear;
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontSize',16);

%% parameter
% Universal constants
F_RT  = 38.68;         % Constant: F/RT, unit 1/V
VLi_F = 1.347*10^-4;   % MolarVolume_Li/F, in unit cm^3/(s*A)
ifg   = 0;             % Figure plot flag
% Electrochemistry parameters
sgmSE = 0.1;       % The ionic conductivity of the SE, mS/cm
% Geometric parameters
hSE   = 20.0;  % thickness of the SE, in unit um
wSE   = 8.0;   % width of the SE, in unit um
rp    = 0.05;  % the round radius at corner point, in unit um
pd    = 0.5;     % pore diameter, in unit um
% Temporal parameters
dt    = 0.1;       % Timestep, in unit hours
nT    = 20;       % Total number of steps
% Variables
hLi   = zeros(nT,1); % The length of Li due to Li growth at each timestep
hLi(1)= 0.2;         % Initial length of Li, in unit um
dltL  = 0.2;         % The segment at Li tip used to compute length growth of each step
nSmp  = 10000;        % Number of sampling points for recoording Li thickness
ySE   = linspace(0,hSE-1.5*rp,nSmp)'; % y Coordinate of sampling points, ending at 1.5*rp below top surface
tLi   = zeros(nSmp, nT);           % The Li thickness at each timestep on sampling points

%% Calculation
% Step 1: Get the intinal Li thickness (dtLi) and length (dL) increasement
[dtLi,dhLi]  = thk_dL(1,dt,wSE,ySE,dltL);

for iT = 2 : nT
% Step 2: compute the current Li thickness and length
    tLi(:,iT) = tLi(:,iT-1) + dtLi;
    hLi(iT)   = hLi(iT-1) + dhLi;

% Step 3: Modify the Li thickness and length constraint by pore size
    [tLi(:,iT), hLi(iT)] = pore_constrain(ySE, pd, tLi(:,iT), hLi(iT));

% Step 4: Generate the new bash script with new Li length
    fid  = fopen('matR.sh','r');
    f    = fread(fid,'*char')';
    fclose(fid);
    f    = strrep(f,'howardtu',num2str(hLi(iT)));
    f    = strrep(f,'tttt',['t',num2str(iT)]);
    fid  = fopen('matRun.sh','w');
    fprintf(fid,'%s',f);
    fclose(fid);
 
% Step 4: run the simulation
    system('bash matRun.sh');

% Step 5: compute the Li thickness (dtLi) and length (dL) increasement for next step
    [dtLi,dhLi]  = thk_dL(iT,dt,wSE,ySE,dltL);
    
end

%save
save tLi;
save hLi;
