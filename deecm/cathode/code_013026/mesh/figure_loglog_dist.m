clear all
clc
close all

load('nmc_constant.mat')

figure('Position',[100 100 600 500]);
cmap = colormap(jet(6));
colormap(cmap)
plot(lps1.dia, 100*lps1.dVol,'-*','Color',cmap(1,:));
hold on
plot(lps2.dia, 100*lps2.dVol,'-*','Color',cmap(2,:));
hold on
plot(lps4.dia, 100*lps4.dVol,'-*','Color',cmap(3,:));
hold on
plot(lps6.dia, 100*lps6.dVol,'-*','Color',cmap(4,:));
hold on
plot(lps8.dia, 100*lps8.dVol,'-*','Color',cmap(5,:));
hold on
plot(lps10.dia, 100*lps10.dVol,'-*','Color',cmap(6,:));
hold on
plot(ncm1.dia, 100*ncm1.dVol,'--ok');
set(gca,'FontSize',20);
xlabel('Diameter ($\mu$m)','Interpreter','latex');
ylabel('Probability','Interpreter','latex');
axis([1 6 0.01 40])

clear all

load('lps_constant.mat')

figure('Position',[100 100 600 500]);
cmap = colormap(jet(6));
colormap(cmap)
plot(ncm1.dia, 100*ncm1.dVol,'--o','Color',cmap(1,:));
hold on
plot(ncm2.dia, 100*ncm2.dVol,'--o','Color',cmap(2,:));
hold on
plot(ncm4.dia, 100*ncm4.dVol,'--o','Color',cmap(3,:));
hold on
plot(ncm6.dia, 100*ncm6.dVol,'--o','Color',cmap(4,:));
hold on
plot(ncm8.dia, 100*ncm8.dVol,'--o','Color',cmap(5,:));
hold on
plot(ncm10.dia, 100*ncm10.dVol,'--o','Color',cmap(6,:));
hold on
plot(lps1.dia, 100*lps1.dVol,'-*k');
set(gca,'FontSize',20);
xlabel('Diameter ($\mu$m)','Interpreter','latex');
ylabel('Probability','Interpreter','latex');
axis([1 6 0.01 40])