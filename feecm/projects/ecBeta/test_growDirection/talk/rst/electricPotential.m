clc; clear; myFigureSetting;

ifg = 0;
xSE   = 500;      % length of the SE, in unit um
ySE   = 500;      % height of the SE, in unit um

SEPot0    = csvread('mdl0_SEInt_potential_0001.csv',1,0);
LiPot0    = csvread('mdl0_LiInt_potential_0001.csv',1,0);


SEPot1    = csvread('mdl1_SEInt_potential_0001.csv',1,0);
LiPot1    = csvread('mdl1_LiInt_potential_0001.csv',1,0);

SEPot2    = csvread('mdl2_SEInt_potential_0001.csv',1,0);
LiPot2    = csvread('mdl2_LiInt_potential_0001.csv',1,0);


% Plot the potential of Li
ifg = ifg + 1;
figure(ifg)
plot(LiPot0(:,1),LiPot0(:,5)*10^6,'-r',LiPot1(:,1),LiPot1(:,5)*10^6,'-b',LiPot2(:,1),LiPot2(:,5)*10^6,'-k')
legend('H_{Li}=0.01um','H_{Li}=0.1um','H_{Li}=1um');

ifg = ifg + 1;
figure(ifg)
plot(SEPot0(:,1),SEPot0(:,5),'-r',SEPot1(:,1),SEPot1(:,5),'-b',SEPot2(:,1),SEPot2(:,5),'-k')
legend('H_{Li}=0.01um','H_{Li}=0.1um','H_{Li}=1um');