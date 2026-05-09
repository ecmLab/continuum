clc; clear;
set(0, 'DefaultLineLineWidth', 2);
set(0,'defaultAxesFontSize',16);

%% parameter
% Universal constants
F_RT  = 38.68;         % Constant: F/RT, unit 1/V
VLi_F = 1.347*10^-4;   % MolarVolume_Li/F, in unit cm^3/(s*A)
ifg   = 0;             % Figure plot flag
% Electrochemistry parameters
sgmSE = 0.1;           % The ionic conductivity of the SE, mS/cm
% Geometric parameters
hSE   = 16.0;  % thickness of the SE, in unit um
wSE   = 10.0;  % radius of the SE, in unit um
dw    = 2;      % width of the crack, in unit um
hLi   = 2;      % length of Li, in unit um
hAg   = hLi/2;  % center location of the Ag, in unit um
dAg   = 0.1;      % radius of the Ag particle, in unit um

%% Load CSV data from MOOSE calculation: structure of the csv data
% load current
% Column 1-2: current density along the x and y direction, in unit mA/cm^2
% Column 4-5: coordinate in the x and y direction, in unit um
% Load the pure SE results with Li metal length of 1um
% Ag1Ref      = csvread(['rst/1AgRef_anode_current_0001.csv'],1,0);
% Ag1Ref      = sortrows(Ag1Ref,[5,4],{'ascend' 'ascend'});       % sort the coordinate with Y and x ascend
% load potential
% Column 2: potential along the li metal, in unit mA/cm^2
% Column 3-4: coordinate in the x and y direction, in unit um
Ag1Ref      = csvread(['rst/1AgRef_anode_potential_0001.csv'],1,0);
Ag1Ref      = sortrows(Ag1Ref,[4,3],{'ascend' 'descend'});       % sort the coordinate with Y and x ascend
% Ag1Ref(Ag1Ref(:,5)==0,:)   = sortrows(Ag1Ref(Ag1Ref(:,5)==0,:),4,'descend');       % sort the coordinate with Y and x ascend
% crnt     = Ag1Ref(:,1:2); % Sort the current density accordingly
crnt     = Ag1Ref(:,1:2); % Sort the current density accordingly
cord        = Ag1Ref(:,4:5); % coordinates of points
% Compute the arch length along Li metal, and the normal current density
archL    = zeros(length(cord),1);
nrCnt    = zeros(length(cord),1);
nrCnt(1) = -crnt(1,2);
for ip = 2 : length(cord)
    archL(ip) = archL(ip-1) + norm(cord(ip,:)-cord(ip-1,:));
    if cord(ip,2) == 0
        nrCnt(ip) = -crnt(ip,2);
%     elseif cord(ip,1) == dw/2
    else
        nrCnt(ip) = -crnt(ip,1);
%     else
%         dNr       = [dw/2,hAg] - cord(ip,:);
%         nrCnt(ip) = tmpCrnt(ip,:)*dNr'/norm(dNr);
    end
end
Ag1Ref = [archL, nrCnt];

% Load the Ag-SE results with Li metal length of 1um and 1 Ag particle with size 0.1um
Ag1  = csvread(['rst/1Ag_anode_current_0001.csv'],1,0);
[cord,indx] = sortrows(Ag1(:,4:5),[1,2],{'descend' 'ascend'}); % sort the coordinate with x ascend and y dscend
crnt        = Ag1(indx,1:2); % Sort the current density accordingly
% Compute the arch length along Li metal, and the normal current density
archL    = zeros(length(cord),1);
nrCnt    = zeros(length(cord),1);
nrCnt(1) = crnt(1,2);
for ip = 2 : length(cord)
    archL(ip) = archL(ip-1) + norm(cord(ip,:)-cord(ip-1,:));
    if cord(ip,2) == 0
        nrCnt(ip) = crnt(ip,2);
    elseif cord(ip,1) == dw/2
        nrCnt(ip) = crnt(ip,1);
    else
        dNr       = cord(ip,:) - [dw/2,hAg];
        nrCnt(ip) = crnt(ip,:)*dNr'/norm(dNr);
    end
end
Ag1 = [archL, nrCnt];

%% Data analysis
ifg = ifg + 1;
figure(ifg)
plot(Ag1Ref(1:end-5,1), abs(Ag1Ref(1:end-5,2)),'--k',Ag1(1:end-5,1), abs(Ag1(1:end-5,2)),'-b')
% plot(Ag1Ref(:,1), Ag1Ref(:,2),'--k')


