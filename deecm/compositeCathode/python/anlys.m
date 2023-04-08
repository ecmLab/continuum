
pi     = 3.1415926;
cn_F   = 96485.3329;
mss_ncm = 96.95406;
Qs_max  = cn_F/mss_ncm * 0.277778;

tmp = load('x_voltage');

mr15_lps3_I05=load('ncm_mr15_lps3_I05.txt');
mr20_lps3_I005=load('ncm_mr20_lps3_I005.txt');
mr20_lps3_I05=load('ncm_mr20_lps3_I05.txt');
mr20_lps5_I05=load('ncm_mr20_lps5_I05.txt');
mr25_lps3_I05=load('ncm_mr25_lps3_I05.txt');
liquidCell  = [tmp(:,1)*Qs_max,tmp(:,2)];

% Interpolate data
qx    = linspace(57.5,176.7,500);
I05   = mr20_lps3_I05(223:end,:);
I005  = mr20_lps3_I005;
lCell = liquidCell(514:end,:);

lclv   = interp1(lCell(:,1),lCell(:,2),qx);
% I05t   = interp1(I05(:,2),I05(:,1),qx)/3600;
I05t   = linspace(0,2.1,500);
I05v   = interp1(I05(:,2),I05(:,3),qx);
% I005t  = interp1(I005(:,2),I005(:,1),qx)/3600;
I005t   = linspace(0,22,500);
I005v  = interp1(I005(:,2),I005(:,3),qx);

figure(1)
plot(qx,lclv,qx,I005v,qx,I05v)
legend('Liquid Cell','I=0.05mA/cm^2','I=0.5mA/cm^2')
xlabel('Capacity (mAh/g)');
ylabel('Voltage (V)');

figure(3)
plot(I005t,(I005v-lclv)*1000,I05t,(I05v-lclv)*1000)
legend('I=0.05mA/cm^2','I=0.5mA/cm^2')
xlabel('Time');
ylabel('Voltage (mV)');

% figure(1)
% plot(liquidCell(:,1),liquidCell(:,2),mr20_lps3_I005(:,2),mr20_lps3_I005(:,3),mr20_lps3_I05(:,2),mr20_lps3_I05(:,3))
% legend('Liquid Cell','I=0.05mA/cm^2','I=0.5mA/cm^2')
% xlabel('Capacity (mAh/g)');
% ylabel('Voltage (V)');
% 
% figure(2)
% plot(liquidCell(:,1),liquidCell(:,2),mr15_lps3_I05(:,2),mr15_lps3_I05(:,3),mr20_lps3_I05(:,2),mr20_lps3_I05(:,3),mr25_lps3_I05(:,2),mr25_lps3_I05(:,3))
% legend('Liquid Cell','Weight ratio 0.7','Weight ratio 0.75','Weight ratio 0.8')
% xlabel('Capacity (mAh/g)');
% ylabel('Voltage (V)');
% 
% figure(3)
% plot(liquidCell(:,1),liquidCell(:,2),mr20_lps3_I05(:,2),mr20_lps3_I05(:,3),mr20_lps5_I05(:,2),mr20_lps5_I05(:,3))
% legend('Liquid Cell','LPS-3um','LPS-5um')
% xlabel('Capacity (mAh/g)');
% ylabel('Voltage (V)');

