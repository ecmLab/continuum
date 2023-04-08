clc;clear;
% load data
lmd   = [0.3333,0.4167,0.625,1,1.667,3.333,4,4.167,6.25,8,12.5];
wt    = [56:1:90]';
aCAM  = load('active_CAM.txt');

% 2D interpolate coordinates
ntp     = 500; % Number of interpolate point
lmd_int = linspace(0.3333,12.5,ntp);
wt_int  = linspace(56,90,ntp)';        % Weight ratio

% Mesh grid all cordinates
lmd_mesh = repmat(lmd,length(wt),1);
wt_mesh  = repmat(wt,1,length(lmd));
lmd_int_mesh = repmat(lmd_int,ntp,1);
wt_int_mesh  = repmat(wt_int,1,ntp);

% Interpolate
%aCAM_int   = interp2(lmd_mesh,wt_mesh,aCAM,lmd_int_mesh,wt_int_mesh,'spline');
aCAM_int   = interp2(lmd_mesh,wt_mesh,aCAM,lmd_int_mesh,wt_int_mesh,'makima');

% Find the lmd-wt curve with the given target active CAM
% tgt = [0.75,0.9, 0.98];
tgt = 0.97;
lmd_wt = zeros(ntp,2);
lmd_wt(:,1) = wt_int;

for i = 1 : length(tgt)
    for j = 1 : ntp
        [tmp,idx] = min(abs(aCAM_int(j,:) - tgt(i)));
        if tmp < 0.003
            lmd_wt(j,i+1) = lmd_int(idx);
        end
    end
end

% Plot
% Initial data
figure(1)
contourf(lmd_int_mesh,wt_int_mesh,aCAM_int,'LineStyle','none','LevelStep',0.001);
%contourf(lmd_mesh,wt_mesh,aCAM,'LineStyle','none','LevelStep',0.001);

% create a default color map ranging from red to light pink
len1 = 400;
len2 = 200;
len3 = 200;
len4 = 75;
cprp = [1, 0, 1]/1.1;
cred = [1, 0, 0];
corg = [1,0.5,0];
cylw = [1, 1, 0]/1.1;
cgrn = [0,150,0]/255;
c_p1 = [linspace(cprp(1),cred(1),len1)', linspace(cprp(2),cred(2),len1)', linspace(cprp(3),cred(3),len1)'];

c_p2 = [linspace(cred(1),corg(1),len2)', linspace(cred(2),corg(2),len2)', linspace(cred(3),corg(3),len2)'];
c_p3 = [linspace(corg(1),cylw(1),len3)', linspace(corg(2),cylw(2),len3)', linspace(corg(3),cylw(3),len3)'];

c_p4 = [linspace(cylw(1),cgrn(1),len4)', linspace(cylw(2),cgrn(2),len4)', linspace(cylw(3),cgrn(3),len4)'];
c_pm = [c_p1;c_p2;c_p3;c_p4];
colormap(c_pm)

% % Interpolated data
% figure(2)
% contourf(lmd_int_mesh,wt_int_mesh,aCAM_int,'LineStyle','none','LevelStep',0.01);
% map = [0.2 0.1 0.5
%     0.15 0.3 0.65
%     0.1 0.5 0.8
%     0.2 0.7 0.6
%     0.8 0.7 0.3
%     0.9 1 0];
% colormap(map)

% Target data
% figure(2)
% idx_plt = lmd_wt(:,2)>0;
% xx = lmd_wt(idx_plt,1);
% yy = lmd_wt(idx_plt,2);
% plot(xx,yy,'-*');
% axis([60 82 0.7 8])

% figure(3)
% for i = 1:length(tgt)
%     idx_plt = lmd_wt(:,i+1)>0;
%     xx = lmd_wt(idx_plt,1);
%     yy = lmd_wt(idx_plt,i+1);
%     if i == 1
%         fp = polyfit(xx,yy,10);
%     elseif i == 2
%         fp = polyfit(xx,yy,11);
%     elseif i == 3
%         fp = polyfit(xx,yy,7);
%     end
%     yp = polyval(fp,xx);
%     plot(xx,yy,'*',xx,yp);
%     hold on
% end
% hold off
% axis([60 82 0.7 8])