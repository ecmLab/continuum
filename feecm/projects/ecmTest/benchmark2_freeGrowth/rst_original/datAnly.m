clc;clear;

% parameter
mdName = 'test';
ifg  = 0;

% Load normal directions of interface boundary
str = strcat('nodal_disp = csvread(''',mdName,'_nodal_disp_0001.csv'',1,0);');
eval(str);
str = strcat('qp_nrml = load(''',mdName,'_swell_normal.txt'');');
eval(str);

%% Output data
% dlmwrite('bndNd_cord.csv',metal_nrml(:,1:3), 'precision',12);

%% Data plot
% plot boundary points with metal normals
% ifg = ifg + 1;
% figure(ifg)
% plot(metal_nrml(:,1),metal_nrml(:,2),'.');
% hold on
% quiver(metal_nrml(:,1),metal_nrml(:,2),metal_nrml(:,5),metal_nrml(:,6),0.2)
% hold off
% axis([-0.01,0.26,-0.21,0.01])

% plot boundary points with normals at ceramic side
ifg = ifg + 1;
figure(ifg)
plot(nodal_disp(:,1),nodal_disp(:,2),'square');
hold on
quiver(qp_nrml(:,1),qp_nrml(:,2),qp_nrml(:,3),qp_nrml(:,4),0.2)
hold off
% axis([-0.02+min(nodal_disp(:,1)),0.02+max(nodal_disp(:,1)),-0.02+min(nodal_disp(:,2)), 0.02+min(nodal_disp(:,1))])
% axis([-0.02,0.52,-0.12,0.02]);
axis square
