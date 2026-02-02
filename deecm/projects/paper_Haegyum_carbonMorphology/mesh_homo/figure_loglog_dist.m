clear all
clc
close all

S = renameLegacyFields(load('nmc_constant.mat'));

figure('Position',[100 100 600 500]);
cmap = colormap(jet(6));
colormap(cmap)
plot(S.carbon1.dia, 100*S.carbon1.dVol,'-*','Color',cmap(1,:));
hold on
plot(S.carbon2.dia, 100*S.carbon2.dVol,'-*','Color',cmap(2,:));
hold on
plot(S.carbon4.dia, 100*S.carbon4.dVol,'-*','Color',cmap(3,:));
hold on
plot(S.carbon6.dia, 100*S.carbon6.dVol,'-*','Color',cmap(4,:));
hold on
plot(S.carbon8.dia, 100*S.carbon8.dVol,'-*','Color',cmap(5,:));
hold on
plot(S.carbon10.dia, 100*S.carbon10.dVol,'-*','Color',cmap(6,:));
hold on
plot(S.drx1.dia, 100*S.drx1.dVol,'--ok');
set(gca,'FontSize',20);
xlabel('Diameter ($\mu$m)','Interpreter','latex');
ylabel('Probability','Interpreter','latex');
axis([1 6 0.01 40])

clear S

S = renameLegacyFields(load('lps_constant.mat'));

figure('Position',[100 100 600 500]);
cmap = colormap(jet(6));
colormap(cmap)
plot(S.drx1.dia, 100*S.drx1.dVol,'--o','Color',cmap(1,:));
hold on
plot(S.drx2.dia, 100*S.drx2.dVol,'--o','Color',cmap(2,:));
hold on
plot(S.drx4.dia, 100*S.drx4.dVol,'--o','Color',cmap(3,:));
hold on
plot(S.drx6.dia, 100*S.drx6.dVol,'--o','Color',cmap(4,:));
hold on
plot(S.drx8.dia, 100*S.drx8.dVol,'--o','Color',cmap(5,:));
hold on
plot(S.drx10.dia, 100*S.drx10.dVol,'--o','Color',cmap(6,:));
hold on
plot(S.carbon1.dia, 100*S.carbon1.dVol,'-*k');
set(gca,'FontSize',20);
xlabel('Diameter ($\mu$m)','Interpreter','latex');
ylabel('Probability','Interpreter','latex');
axis([1 6 0.01 40])

function out = renameLegacyFields(data)
    out = data;
    fields = fieldnames(data);
    for i = 1:numel(fields)
        f = fields{i};
        newName = f;
        if startsWith(f,'ncm')
            newName = strrep(f,'ncm','drx');
        elseif startsWith(f,'lps')
            newName = strrep(f,'lps','carbon');
        end
        if ~strcmp(newName,f)
            out.(newName) = data.(f);
        end
    end
end
