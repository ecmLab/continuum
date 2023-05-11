clc; clear; myFigureSetting;

ifg = 0;
xSE   = 20;      % width of the SE, in unit um
ySE   = 40;      % thickness of the SE, in unit um
dw   = 2;        % width of the defect located at the top middle, in unit um
dh   = 2;        % length of the defect, in unit um

%Load anode potential
andPot_sgm001    = csvread('../rst/andPot_sgm-2.csv',1,0);
andPot_sgm01     = csvread('../rst/andPot_sgm-1.csv',1,0);
andPot_sgm0      = csvread('../rst/andPot_sgm0.csv',1,0);
andPot_sgm1      = csvread('../rst/andPot_sgm1.csv',1,0);
%Load anode current
andCrnt_sgm001    = csvread('../rst/andCrnt_sgm-2.csv',1,0);
andCrnt_sgm01     = csvread('../rst/andCrnt_sgm-1.csv',1,0);
andCrnt_sgm0      = csvread('../rst/andCrnt_sgm0.csv',1,0);
andCrnt_sgm1      = csvread('../rst/andCrnt_sgm1.csv',1,0);

% Plot the potential
ifg = ifg + 1;
figure(ifg)
plot(andPot_sgm001(:,3),andPot_sgm001(:,2),'-r',andPot_sgm01(:,3),andPot_sgm01(:,2),'-b');
hold on
plot(andPot_sgm0(:,3),andPot_sgm0(:,2),'-k',andPot_sgm1(:,3),andPot_sgm1(:,2),'-g');
hold off
legend('0.01mS/cm','0.1mS/cm','1mS/cm','10mS/cm');
