clc; clear; myFigureSetting;

ifg = 0;
xSE   = 500;      % length of the SE, in unit um
ySE   = 500;      % height of the SE, in unit um

depRt0 = load('mdl0.txt');
depRt1 = load('mdl1.txt');
depRt2 = load('mdl2.txt');
depRt4 = load('mdl4.txt');

% Plot the potential of Li
ifg = ifg + 1;
figure(ifg)
plot(depRt0(:,1),depRt0(:,2),'-r',depRt1(:,1),depRt1(:,2),'-b',depRt2(:,1),depRt2(:,2),'-k')
legend('Li-0.01um','Li-0.1um','Li-1um');

ifg = ifg + 1;
figure(ifg)
plot(depRt2(:,1)+500,depRt2(:,2),'-b',depRt4(:,1)+500,depRt4(:,2),'-k')
legend('L_{Li}=5um','L_{Li}=40um');