clc; clear; myFigureSetting;

ifg = 0;
xSE   = 500;      % length of the SE, in unit um
ySE   = 500;      % height of the SE, in unit um

SEInt0_pot    = csvread('mdl0_SEInt_potential_0001.csv',1,0);
tmp          = sortrows(SEInt0_pot(SEInt0_pot(:,1)==0,:),2);
SEInt0_pot(1:length(tmp),:) = tmp;
LiInt0_pot    = csvread('mdl0_LiInt_potential_0001.csv',1,0);
tmp          = sortrows(LiInt0_pot(LiInt0_pot(:,1)==0,:),2);
LiInt0_pot(1:length(tmp),:) = tmp;

SEInt1_pot    = csvread('mdl1_SEInt_potential_0001.csv',1,0);
tmp          = sortrows(SEInt1_pot(SEInt1_pot(:,1)==0,:),2);
SEInt1_pot(1:length(tmp),:) = tmp;
LiInt1_pot    = csvread('mdl1_LiInt_potential_0001.csv',1,0);
tmp          = sortrows(LiInt1_pot(LiInt1_pot(:,1)==0,:),2);
LiInt1_pot(1:length(tmp),:) = tmp;

SEInt2_pot    = csvread('mdl2_SEInt_potential_0001.csv',1,0);
tmp          = sortrows(SEInt2_pot(SEInt2_pot(:,1)==0,:),2);
SEInt2_pot(1:length(tmp),:) = tmp;
LiInt2_pot    = csvread('mdl2_LiInt_potential_0001.csv',1,0);
tmp          = sortrows(LiInt2_pot(LiInt2_pot(:,1)==0,:),2);
LiInt2_pot(1:length(tmp),:) = tmp;


% Plot the potential of Li
ifg = ifg + 1;
figure(ifg)
LiPot0 = LiInt0_pot(LiInt0_pot(:,2)==ySE,:);
LiPot1 = LiInt1_pot(LiInt1_pot(:,2)==ySE,:);
LiPot2 = LiInt2_pot(LiInt2_pot(:,2)==ySE,:);

plot(LiPot0(:,1),LiPot0(:,5),'-or',LiPot1(:,1),LiPot1(:,5),'-b',LiPot2(:,1),LiPot2(:,5),'-k')
legend('Li-0.01um','Li-0.5um','Li-2um');

ifg = ifg + 1;
figure(ifg)
SEPot0 = SEInt0_pot(SEInt0_pot(:,2)==ySE,:);
SEPot1 = SEInt1_pot(SEInt1_pot(:,2)==ySE,:);
SEPot2 = SEInt2_pot(SEInt2_pot(:,2)==ySE,:);

plot(SEPot0(:,1),SEPot0(:,5),'-or',SEPot1(:,1),SEPot1(:,5),'-b',SEPot2(:,1),SEPot2(:,5),'-k')
legend('Li-0.01um','Li-0.5um','Li-2um');