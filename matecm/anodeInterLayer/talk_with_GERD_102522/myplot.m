tmp = load('data/vol_20.mat');
vol20  = tmp.vol_20;
tmp = load('data/vol_40.mat');
vol40  = tmp.vol_40;
tmp = load('data/vol_80.mat');
vol80  = tmp.vol_80;

tmp = load('data/narrow_11.mat');
narrow11  = tmp.aa;
tmp = load('data/narrow_30.mat');
narrow30  = tmp.aa;
tmp = load('data/wide_11.mat');
wide11  = tmp.aa;
tmp = load('data/wide_30.mat');
wide30  = tmp.aa;

stableEqm = 3.7622;

figure(1)
plot(vol20(:,1)*100,vol20(:,2),'-*b',vol40(:,1)*100,vol40(:,2),'-xr',vol80(:,1)*100,vol80(:,2),'-og')
yline(stableEqm,'--k');
xlabel('Vol% of Ag');
ylabel('stable IL voltage (mV)');
legend('20nm','40nm','80nm');

figure(2)
plot(vol20(:,1)*100,vol20(:,3),'-*b',vol40(:,1)*100,vol40(:,3),'-xr',vol80(:,1)*100,vol80(:,3),'-og')
xlabel('Vol% of Ag');
ylabel('stable time (seconds)');
legend('20nm','40nm','80nm');

figure(3)
plot(narrow11(:,1),narrow11(:,2),'-b',narrow30(:,1),narrow30(:,2),'-r',wide11(:,1),wide11(:,2),'-g',wide30(:,1),wide30(:,2),'-k')
xlabel('Charging time (Hours)');
ylabel('Voltage (mV)');
legend('Narrow-1.1v%','Narrow-3v%','Wide-1.1v%','Wide-3v%');
axis([0,0.1,0,100]);