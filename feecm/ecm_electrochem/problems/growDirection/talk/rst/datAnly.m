clc; clear; myFigureSetting;

%% parameter
% Geometry of the model
str = 'mdl4';
xSE   = 500;      % length of the SE, in unit um
ySE   = 500;      % height of the SE, in unit um
rp    = 0.1;   % The round radius
xLi   = 5;    % length of the Li, in unit um
yLi   = 1;  % height of the Li, in unit um

if strcmp(str,'mdl4')
    xLi   = 40;    % length of the Li, in unit um
elseif strcmp(str,'mdl5')
    xLi   = 0.1;   % length of the Li, in unit um
end
if strcmp(str,'mdl0')
    yLi   = 0.01;   % height of the Li, in unit um
elseif strcmp(str,'mdl1')
    yLi   = 0.1;   % height of the Li, in unit um
    
end

% Electrochemistry parameters
i_exc = 1.3;       % The exchange current density of lps, mA/cm^2
alfa  = 0.5;       % Reaction rate
i_applied = 1.0;  % The applied current density
F_RT  = 38.68; % Constant: F/RT, unit 1/V
ifg   = 0;      % Figure plot flag
VLi_F = 1.3;      % V_Li/F, in unit nm/s

%% Load data
% Total current at each boundary, in unit nA
% Column 1: Total current flow through Li_anode;   Column 2: Total current flow through SE_anode;
% Column 3: Total current flow through SE_cathode; Column 4: Total current flow through interface;
% positive means flow out from the block, negative means flow into the block.
Tot_crnt     = csvread([str,'.csv'],1,1);
% SE potential and current
SEAnode_pot  = csvread([str,'_SEAnode_potential_0001.csv'],1,0);
SEInt_pot    = csvread([str,'_SEInt_potential_0001.csv'],1,0);
SETop_pot    = csvread([str,'_SETop_potential_0001.csv'],1,0);
SECathode_pot= csvread([str,'_SECathode_potential_0001.csv'],1,0);
SEBtm_pot    = csvread([str,'_SEBottom_potential_0001.csv'],1,0);
SEAnode_crnt = csvread([str,'_SEAnode_current_0001.csv'],1,0);
SEAnode_crnt = sortrows(SEAnode_crnt,2);
SEInt_crnt   = csvread([str,'_SEInt_current_0001.csv'],1,0);
% Li potential and current
LiCathode_pot = csvread([str,'_LiCathode_potential_0001.csv'],1,0);
LiInt_pot     = csvread([str,'_LiInt_potential_0001.csv'],1,0);
LiTop_pot     = csvread([str,'_LiTop_potential_0001.csv'],1,0);
LiAnode_crnt  = csvread([str,'_LiAnode_current_0001.csv'],1,0);
LiInt_crnt    = csvread([str,'_LiInt_current_0001.csv'],1,0);

%% Data process
% SE Potential distribution along the edge: SE_left -- Interface -- SE_top -- SE_right -- SE_bottom
nEdg = length(SEAnode_pot)+length(SEInt_pot)+length(SETop_pot)+length(SECathode_pot)+length(SEBtm_pot); %Total number of edge points
SEEdge_pot = zeros(nEdg,4);
SEEdge_pot(:,1:2) = [SEAnode_pot(:,1:2);SEInt_pot(:,1:2);SETop_pot(:,1:2);flipud(SECathode_pot(:,1:2));flipud(SEBtm_pot(:,1:2))];
tmp = vecnorm(SEEdge_pot(2:nEdg,1:2) - SEEdge_pot(1:nEdg-1,1:2),2,2);
% Arch length
SEEdge_pot(1,3) = SEEdge_pot(1,2);
for i = 2 : nEdg
    SEEdge_pot(i,3) = SEEdge_pot(i-1,3) + tmp(i-1);
end
SEEdge_pot(:,4) = [SEAnode_pot(:,5);SEInt_pot(:,5);SETop_pot(:,5);flipud(SECathode_pot(:,5));flipud(SEBtm_pot(:,5))];

% Analytic solution of deposition rate along edge "SE_left -- Interface" according B-V equation
eta = [SEAnode_pot(:,5);SEInt_pot(:,5) - LiInt_pot(:,5)];  % Overpotential, the Li metal voltage is very close to zero (10^-10)
depRtAly = i_exc * (exp(alfa*F_RT*eta) - exp(-(1-alfa)*F_RT*eta));
% Numerical solution of deposition rate along edge "SE_left -- Interface"
depRtNrc = -SEAnode_crnt(:,5);
for i = 1 : length(SEInt_pot)
    if SEInt_crnt(i,1) < rp
        depRtNrc(length(SEAnode_pot)+i) = SEInt_crnt(i,5)*(SEInt_crnt(i,1)-rp)/rp + SEInt_crnt(i,6)*(SEInt_crnt(i,2)-ySE+rp)/rp;
    else
        depRtNrc(length(SEAnode_pot)+i) = SEInt_crnt(i,6);
    end
end

% Numerical solution of deposition rate at interface for SE and Li metal
depRtInt      = zeros(length(SEInt_pot),3);
depRtInt(:,1) = SEEdge_pot(length(SEAnode_pot)+1:length(depRtNrc),3);
depRtInt(:,2) = depRtNrc(length(SEAnode_pot)+1:length(depRtNrc));
for i = 1 : length(SEInt_pot)
    if LiInt_crnt(i,1) < rp
        depRtInt(i,3) = -LiInt_crnt(i,5)*(LiInt_crnt(i,1)-rp)/rp - LiInt_crnt(i,6)*(LiInt_crnt(i,2)-ySE+rp)/rp;
    else
        depRtInt(i,3) = -LiInt_crnt(i,6);
    end
end

% Percentage of electrons flow into the Li metal block to the total applied current
Tot_crnt(1)/Tot_crnt(3)

%% Plot
% Plot SE potential along thickness
% ifg = ifg + 1;
% figure(ifg)
% plot(SEInt_pot(:,1),SEInt_pot(:,5),'.b',SETop_pot(:,1),SETop_pot(:,5),'-k')
% legend('SE-pot-Top','SE-pot-Bottom');
% % Plot the potential of SE at the anode
% ifg = ifg + 1;
% figure(ifg)
% plot(SEAnode_pot(:,2),SEAnode_pot(:,5),'-k')
% legend('SE-pot-Anode');
% % Plot the normal current of SE at the anode
% ifg = ifg + 1;
% figure(ifg)
% plot(SEAnode_crnt(:,2),-SEAnode_crnt(:,5),'-k')
% legend('SE-crnt-Anode');

% Plot the potential of Li along thickness
% ifg = ifg + 1;
% figure(ifg)
% plot(LiInt_pot(:,1),LiInt_pot(:,5),'o',LiTop_pot(:,1),LiTop_pot(:,5),'-k')
% legend('Li-pot-Bottom','Li-pot-Top');
% % Plot the potential of Li at the cathode
% ifg = ifg + 1;
% figure(ifg)
% plot(LiCathode_pot(:,2),LiCathode_pot(:,5),'-k')
% legend('Li-pot-Cathode');
% % Plot the normal current of Li at the anode
% ifg = ifg + 1;
% figure(ifg)
% plot(LiAnode_crnt(:,2),LiAnode_crnt(:,5),'-k')
% legend('Li-crnt-Anode');

% Plot SE potential along the edge: SE_left -- Interface -- SE_top -- SE_right -- SE_bottom
ifg = ifg + 1;
figure(ifg)
plot(SEEdge_pot(:,3),SEEdge_pot(:,4),'.b')
hold on
line([0,SEEdge_pot(length(SEAnode_pot)+length(SEInt_pot),3)*1.1],[1,1]*SEEdge_pot(length(SEAnode_pot)+length(SEInt_pot),4),'Color','black')
% line([1,1]*SEEdge_pot(length(SEAnode_pot)+length(SEInt_pot),3),[min(SEEdge_pot(:,4)),SEEdge_pot(length(SEAnode_pot)+length(SEInt_pot),4)],'Color','black')
hold off
axis([0,max(SEEdge_pot(:,3)),min(SEEdge_pot(:,4)),max(SEEdge_pot(:,4))])
xlabel('Arch length (\mum)')
ylabel('Li+ potential (V)')

% Plot the normal current of SE along edge "SE_left -- Interface", compare with analytic result
ifg = ifg + 1;
figure(ifg)
% plot(SEEdge_pot(1:length(depRtNrc),3),depRtNrc, SEEdge_pot(1:length(depRtAly),3), depRtAly, '.b')
% legend('numerical','Analytic');
plot(SEEdge_pot(1:length(depRtAly),3), depRtAly, '.b')
xlabel('Arch length (\mum)')
ylabel('Deposition rate along interface (mA/cm^2)')

% save data
depRt = [SEEdge_pot(1:length(depRtAly),3),depRtAly*VLi_F];
fid=fopen([str,'.txt'],'w');
fprintf(fid,'%f %f\n', depRt');
fclose(fid)
% Plot the normal current of SE along edge "SE_left -- Interface", compare with analytic result
% ifg = ifg + 1;
% figure(ifg)
% plot(depRtInt(:,1),depRtInt(:,2),'.k', depRtInt(:,1),depRtInt(:,3), '.b')
% legend('SE','Li');
