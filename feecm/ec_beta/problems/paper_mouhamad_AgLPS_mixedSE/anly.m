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
hSE   = 20.0;  % thickness of the SE, in unit um
wSE   = 8.0;  % radius of the SE, in unit um
rp    = 0.05;  % the round radius at corner point, in unit um
hLi   = 2;      % length of Li, in unit um
hAg   = hLi/2;  % center location of the Ag, in unit um
dAg   = 0.1;      % radius of the Ag particle, in unit um

%% Load CSV data from MOOSE calculation: structure of the csv data
% Column 1-2: current density along the x and y direction, in unit mA/cm^2
% Column 4-5: coordinate in the x and y direction, in unit um
% Load the pure SE results with Li metal length of 1um
Ag0      = csvread(['rst/0Ag_anode_current_0001.csv'],1,0);
Agp      = csvread(['rst/0Ag_anode_potential_0001.csv'],1,0);
Ag0      = sortrows(Ag0,[4,5],{'ascend' 'descend'});       % sort the coordinate with Y and x ascend
Agp      = sortrows(Agp,[3,4],{'ascend' 'descend'});       % sort the coordinate with Y and x ascend
crnt     = Ag0(:,1:2); % Sort the current density accordingly
cord     = Ag0(:,4:5); % coordinates of points
% Compute the arch length along Li metal, and the normal current density
archL    = zeros(length(cord),1);
nrCnt    = zeros(length(cord),1);
nrCnt(1) = crnt(1,2);
for ip = 2 : length(cord)
    archL(ip) = archL(ip-1) + norm(cord(ip,:)-cord(ip-1,:));
    if cord(ip,2) == hSE
        nrCnt(ip) = crnt(ip,2);
    elseif cord(ip,1) == wSE
        nrCnt(ip) = crnt(ip,1);
    else
        dNr       = cord(ip,:) - [wSE-rp,hSE-rp];
        nrCnt(ip) = crnt(ip,:)*dNr'/norm(dNr);
    end
end
Ag0 = [archL, nrCnt];

% % Load the Ag-SE results with Li metal length of 1um and 1 Ag particle with size 0.1um
% Ag1  = csvread(['rst/1Ag_anode_current_0001.csv'],1,0);
% [cord,indx] = sortrows(Ag1(:,4:5),[1,2],{'descend' 'ascend'}); % sort the coordinate with x ascend and y dscend
% tmpCrnt        = Ag1(indx,1:2); % Sort the current density accordingly
% % Compute the arch length along Li metal, and the normal current density
% archL    = zeros(length(cord),1);
% nrCnt    = zeros(length(cord),1);
% nrCnt(1) = tmpCrnt(1,2);
% for ip = 2 : length(cord)
%     archL(ip) = archL(ip-1) + norm(cord(ip,:)-cord(ip-1,:));
%     if cord(ip,2) == 0
%         nrCnt(ip) = tmpCrnt(ip,2);
%     elseif cord(ip,1) == dw/2
%         nrCnt(ip) = tmpCrnt(ip,1);
%     else
%         dNr       = cord(ip,:) - [dw/2,hAg];
%         nrCnt(ip) = tmpCrnt(ip,:)*dNr'/norm(dNr);
%     end
% end
% Ag1 = [archL, nrCnt];

%% Data analysis
ifg = ifg + 1;
figure(ifg)
% plot(Ag1Ref(1:end-5,1), abs(Ag1Ref(1:end-5,2)),'--k',Ag1(1:end-5,1), abs(Ag1(1:end-5,2)),'-b')
% plot(Ag0(:,1), Ag0(:,2),'--k')
plot(Ag0(:,1), Agp(:,2),'--k')

