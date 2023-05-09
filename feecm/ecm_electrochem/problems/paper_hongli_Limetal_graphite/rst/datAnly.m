clc; clear; myFigureSetting;

%% parameter
% Geometry of the model
xGr   = 26.55;      % thickness of the GR, in unit um
yGr   = 10;      % height of the Gr, in unit um
xLi   = 5;    % length of the Li, in unit um
vLi   = 0.1;     % The extrusion rate of Li, in unit um/s
% Electrochemistry parameters
i_exc = 1.3;       % The exchange current density of lps, mA/cm^2
alfa  = 0.5;       % Reaction rate
i_applied = 1.0;  % The applied current density
F_RT  = 38.68; % Constant: F/RT, unit 1/V
VLi_F = 1.3;      % V_Li/F, in unit nm/s
En    = 20;      % Li nucleation energy, mV
% Simulation parameter
lLi = linspace(4,18,8);   % Tje Li length, in unit um

%% Load data
%1. Li concentration data at the mechanochemistry phase from lLi = 4:18 um
cLi  = zeros(length(lLi),7);
for iLi = 1 : length(lLi)
	icLi = csvread(['t',num2str(lLi(iLi)),'.csv'],1,0);
    cLi(iLi,:) = icLi(end,:);
end
LiC    = cLi(:,end:-1:2)*6;       % Convert from LixC to LixC6
%2. Li concentration data at lLi=2um, corresponding to stack pressure=10
cLi_10 = csvread(['t2.csv'],1,0);
LiC_10 = cLi_10(:,end:-1:2)*6;    % Convert from LixC to LixC6
%3. Li concentration data at lLi=20um, corresponding to stack pressure=300
cLi_300 = csvread(['t20.csv'],1,0);
LiC_300 = cLi_300(:,end:-1:2)*6;   % Convert from LixC to LixC6
%4. Li concentration data at the electrochemistry phase
cLi_ec = csvread(['t20_ec.csv'],1,0);
LiC_ec = cLi_ec(:,end:-1:2)*6;    % Convert from LixC to LixC6
%5. Voltage profile of the Li-Gr of t24
tmp     = csvread(['t24_miec_potLi.csv'],1,0);
vlt_t24 = tmp(:,[3,4,2]);
% sample voltage data with meshgrid
xx = linspace(0,max(vlt_t24(:,1)));
yy = linspace(0,max(vlt_t24(:,2)));
[xsmp, ysmp]  =  meshgrid(xx,yy);
vltFit        = griddata(vlt_t24(:,1),vlt_t24(:,2),vlt_t24(:,3),xsmp,ysmp)*1000;  % Fit voltage data, in unit mV
[C,h]         = contour(xsmp,ysmp,(max(vltFit(:))-vltFit),[1,1]*En);
xyEn          = C(:,2:end);
%30
tmp     = csvread(['t24_miec30_potLi.csv'],1,0);
vlt_t24 = tmp(:,[3,4,2]);
% sample voltage data with meshgrid
xx = linspace(0,max(vlt_t24(:,1)));
yy = linspace(0,max(vlt_t24(:,2)));
[xsmp, ysmp]   =  meshgrid(xx,yy);
vltFit  = griddata(vlt_t24(:,1),vlt_t24(:,2),vlt_t24(:,3),xsmp,ysmp)*1000;  % Fit voltage data, in unit mV
[C,h]         = contour(xsmp,ysmp,(max(vltFit(:))-vltFit),[1,1]*En);
xyEn30        = C(:,2:end);
%60
tmp     = csvread(['t24_miec60_potLi.csv'],1,0);
vlt_t24 = tmp(:,[3,4,2]);
% sample voltage data with meshgrid
xx = linspace(0,max(vlt_t24(:,1)));
yy = linspace(0,max(vlt_t24(:,2)));
[xsmp, ysmp]   =  meshgrid(xx,yy);
vltFit  = griddata(vlt_t24(:,1),vlt_t24(:,2),vlt_t24(:,3),xsmp,ysmp)*1000;  % Fit voltage data, in unit mV
[C,h]         = contour(xsmp,ysmp,(max(vltFit(:))-vltFit),[1,1]*En);
xyEn60        = C(:,2:end);
%90
tmp     = csvread(['t24_miec90_potLi.csv'],1,0);
vlt_t24 = tmp(:,[3,4,2]);
% sample voltage data with meshgrid
xx = linspace(0,max(vlt_t24(:,1)));
yy = linspace(0,max(vlt_t24(:,2)));
[xsmp, ysmp]   =  meshgrid(xx,yy);
vltFit  = griddata(vlt_t24(:,1),vlt_t24(:,2),vlt_t24(:,3),xsmp,ysmp)*1000;  % Fit voltage data, in unit mV
[C,h]         = contour(xsmp,ysmp,(max(vltFit(:))-vltFit),[1,1]*En);
xyEn90        = C(:,2:end);
%120
tmp     = csvread(['t24_miec120_potLi.csv'],1,0);
vlt_t24 = tmp(:,[3,4,2]);
% sample voltage data with meshgrid
xx = linspace(0,max(vlt_t24(:,1)));
yy = linspace(0,max(vlt_t24(:,2)));
[xsmp, ysmp]   =  meshgrid(xx,yy);
vltFit  = griddata(vlt_t24(:,1),vlt_t24(:,2),vlt_t24(:,3),xsmp,ysmp)*1000;  % Fit voltage data, in unit mV
[C,h]         = contour(xsmp,ysmp,(max(vltFit(:))-vltFit),[1,1]*En);
xyEn120        = C(:,2:end);

%6. voltage
tmp     = csvread(['t10_miec_potLi.csv'],1,0);
vlt_t10 = tmp(tmp(:,2)>0,2:3);
tmp     = csvread(['t10_miec_current.csv'],1,0);
crt_t10 = tmp(tmp(:,4)>0,[2,4]);
tmp     = csvread(['t20_miec_potLi.csv'],1,0);
vlt_t20 = tmp(tmp(:,2)>0,2:3);
tmp     = csvread(['t20_miec_current.csv'],1,0);
crt_t20 = tmp(tmp(:,4)>0,[2,4]);
tmp     = csvread(['t22_miec_potLi.csv'],1,0);
vlt_t22 = tmp(tmp(:,2)>0,2:3);
tmp     = csvread(['t22_miec_current.csv'],1,0);
crt_t22 = tmp(tmp(:,4)>0,[2,4]);
tmp     = csvread(['t24_miec_potLi.csv'],1,0);
vlt_t24 = tmp(tmp(:,2)>0,2:3);
tmp     = csvread(['t24_miec_current.csv'],1,0);
crt_t24 = tmp(tmp(:,4)>0,[2,4]);


%% Plot
ifg   = 0;      % Figure plot flag

% plot for the y value at differnt time
ifg = ifg + 1;
figure(ifg)
sid = [1:6];
plot(LiC(5,:),sid,'--b*',LiC(6,:),sid,'--ko',LiC(7,:),sid,'--rx',LiC(8,:),sid,'--g.')

% contour plot for the y value at stack pressure 10
ifg = ifg + 1;
figure(ifg)
tt = cLi_10(:,1)/60;
[ts,scn] = meshgrid(tt,[1:6]');
contourf(ts,scn,LiC_10');

% contour plot for the y value at stack pressure 300
ifg = ifg + 1;
figure(ifg)
tt = cLi_300(:,1)/60;
[ts,scn] = meshgrid(tt,[1:6]');
contourf(ts,scn,LiC_300');

% Plot voltage
ifg = ifg + 1;
figure(ifg)
[C,h]=contourf(xsmp,ysmp,(max(vltFit(:))-vltFit));
clabel(C,h);
ifg = ifg + 1;
figure(ifg)
plot(xyEn(2,:),xyEn(1,:),xyEn30(2,:),xyEn30(1,:),xyEn60(2,:),xyEn60(1,:),xyEn90(2,:),xyEn90(1,:),xyEn120(2,:),xyEn120(1,:))
set(gca, 'YDir','reverse')
legend('0','30','60','90','120');
% plot(vlt_t20(:,1)*1000,vlt_t20(:,2))
% hold on
% plot(vlt_t22(:,1)*1000,vlt_t22(:,2))
% hold on
% plot(vlt_t24(:,1)*1000,vlt_t24(:,2))
% hold off
% set(gca, 'YDir','reverse')

% Plot current
ifg = ifg + 1;
figure(ifg)
plot(crt_t20(:,2),crt_t20(:,1))
hold on
plot(crt_t22(:,2),crt_t22(:,1))
hold on
plot(crt_t24(:,2),crt_t24(:,1))
hold off