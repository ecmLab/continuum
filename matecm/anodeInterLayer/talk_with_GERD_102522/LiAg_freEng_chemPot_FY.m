% clc;clear;
% set(0, 'DefaultLineLineWidth', 2);
% set(0,'defaultAxesFontSize',16);

AA      = [-0.00187215808,  -1.131459,      4.40865164,     -16.7302443,    33.5012958,     -29.9808132,    9.93696756;...
           -0.03250496,     -0.91556805,    2.81440608,     -9.88229831,    20.07645277,    -18.08822859,   6.02541968;...
            430.68417051,   -3275.69649334, 10322.11796142, -17261.4968529, 16158.33241541, -8027.84968135, 1653.91546536];
tang_xp = [0.366,0.405, 0.562,0.637, 0.672,0.766, 0.825,0.884];   % The tangent points of the Gibbs free energy between phases
tang_yp = tang_xp ./ (1-tang_xp);                        % The tangent points of the equilibrium potential between phases
Npp     = 10000;

% The free energy of three phases, in unit J/mol
xp          = linspace(0.0,0.999,Npp)';
freEng_bcc  = freEng_AgLi(xp,AA,1);
freEng_fcc  = freEng_AgLi(xp,AA,2);
freEng_gma  = freEng_AgLi(xp,AA,3);
% The Li chemical potential of each phase, in unit J/mol
chemPot_bcc = chemPot_AgLi(xp,AA,1);
chemPot_fcc = chemPot_AgLi(xp,AA,2);
chemPot_gma = chemPot_AgLi(xp,AA,3);

% The Gibbs free energy and Li chemical potential of the AgLiy with overshoot
freEng_overShoot  = zeros(Npp,1);
chemPot_overShoot = zeros(Npp,1);
for i = 1 : Npp
    if xp(i) < 0.9
        [freEng_overShoot(i),idx] = min([freEng_bcc(i),freEng_fcc(i),freEng_gma(i)]);
        eqm = [chemPot_bcc(i),chemPot_fcc(i),chemPot_gma(i)];
        chemPot_overShoot(i) = eqm(idx);
    else
        chemPot_overShoot(i) = chemPot_bcc(i);
    end
end

% The Gibbs free energy and Li chemical potential of the AgLiy at thermodynamic equilibrium
freEng_thermo  = freEng_overShoot;
chemPot_thermo = chemPot_overShoot;
% Find the co-tangent line between bcc and fcc phase
idx_fcc_bcc                  = find(xp>tang_xp(1) & xp<tang_xp(2));
chemPot_thermo(idx_fcc_bcc)  = chemPot_thermo(idx_fcc_bcc(1)-1);
dedx                         = (freEng_thermo(idx_fcc_bcc(end))-freEng_thermo(idx_fcc_bcc(1)))/(xp(idx_fcc_bcc(end))-xp(idx_fcc_bcc(1)));
freEng_thermo(idx_fcc_bcc)   = dedx * (xp(idx_fcc_bcc)-xp(idx_fcc_bcc(1))) + freEng_thermo(idx_fcc_bcc(1));
% Find the co-tangent line between bcc and gamma1 phase
idx_bcc_gma1                 = find(xp>tang_xp(3) & xp<tang_xp(4));
chemPot_thermo(idx_bcc_gma1) = chemPot_thermo(idx_bcc_gma1(1)-1);
dedx                         = (freEng_thermo(idx_bcc_gma1(end))-freEng_thermo(idx_bcc_gma1(1)))/(xp(idx_bcc_gma1(end))-xp(idx_bcc_gma1(1)));
freEng_thermo(idx_bcc_gma1)  = dedx * (xp(idx_bcc_gma1)-xp(idx_bcc_gma1(1))) + freEng_thermo(idx_bcc_gma1(1));
% Find the co-tangent line between gamma1 and gamma2 phase
idx_gma1_gma2                = find(xp>tang_xp(5) & xp<tang_xp(6));
chemPot_thermo(idx_gma1_gma2)= chemPot_thermo(idx_gma1_gma2(2)-1);
dedx                         = (freEng_thermo(idx_gma1_gma2(end))-freEng_thermo(idx_gma1_gma2(1)))/(xp(idx_gma1_gma2(end))-xp(idx_gma1_gma2(1)));
freEng_thermo(idx_gma1_gma2) = dedx * (xp(idx_gma1_gma2)-xp(idx_gma1_gma2(1))) + freEng_thermo(idx_gma1_gma2(1));
% Find the co-tangent line between gamma2 and bcc Li phase
idx_gma2_bcc                 = find(xp>tang_xp(7) & xp<tang_xp(8));
chemPot_thermo(idx_gma2_bcc) = chemPot_thermo(idx_gma2_bcc(1)-1);
dedx                         = (freEng_thermo(idx_gma2_bcc(end))-freEng_thermo(idx_gma2_bcc(1)))/(xp(idx_gma2_bcc(end))-xp(idx_gma2_bcc(1)));
freEng_thermo(idx_gma2_bcc)  = dedx * (xp(idx_gma2_bcc)-xp(idx_gma2_bcc(1))) + freEng_thermo(idx_gma2_bcc(1));

% The Li chemical potential of the AgLiy when voltage overshoot to zero from BCC phase;
chemPot_bccShoot = chemPot_thermo;
for i = 1 : Npp
    if xp(i) > tang_xp(3)
        chemPot_bccShoot(i) = chemPot_bcc(i);
    end
end
idtp = find(chemPot_bccShoot>0);
chemPot_bccShoot(idtp) = 0;

% The Li chemical potential of the AgLiy when voltage overshoot to zero from gama phase;
chemPot_gmaShoot = chemPot_thermo;
for i = 1 : Npp
    if xp(i) > tang_xp(7)
        chemPot_gmaShoot(i) = chemPot_gma(i);
    end
end
idtp = find(chemPot_gmaShoot>0);
chemPot_gmaShoot(idtp) = 0;

%% Plot
% 1. Plot the free energy of each phase
ifg = ifg + 1;
figure(ifg)
plot(xp,freEng_bcc,'-.',xp,freEng_fcc,'-.',xp(Npp*0.62:end),freEng_gma(Npp*0.62:end),'--');
xlabel('x in Ag_{1-x}Li_x')
ylabel('Gibbs free energy (J/mol)')
legend('BCC','FCC','GAMMA');
% 2. Plot the overall free energy and equilibrium voltage vs. x in Ag_{1-x}Li_x
% ifg = ifg + 1;
% figure(ifg)
% plot(xp,freEng_overShoot,xp,freEng_thermo,'--');
% xlabel('x in Ag_{1-x}Li_x')
% ylabel('Gibbs free energy (J/mol)')
% legend('Overshoot','Thermo');
ifg = ifg + 1;
figure(ifg)
plot(xp,chemPot_overShoot,xp,chemPot_thermo,'--');
xlabel('x in Ag_{1-x}Li_x')
ylabel('Li chemical potential (J/mol)')
legend('Overshoot','Thermo');

% 3. Plot the overall free energy and equilibrium voltage vs. y in AgLi_y
% ifg = ifg + 1;
% figure(ifg)
% ypp = xp ./ (1-xp);
% plot(ypp,freEng_overShoot./(1-xp),ypp,freEng_thermo./(1-xp),'--');
% xlabel('y in AgLi_y')
% ylabel('Gibbs free energy (J/mol)')
% legend('Overshoot','Thermo');
ifg = ifg + 1;
figure(ifg)
ypp = xp ./ (1-xp);
% plot(yp,chemPot_bcc,'--',yp,chemPot_fcc,'--',yp(Np*0.65:end),chemPot_gma(Np*0.65:end),'--',yp,chemPot_overShoot,yp,chemPot_thermo);
plot(xp,-chemPot_overShoot/FF*1000,xp,-chemPot_thermo/FF*1000,'--',xp,-chemPot_bccShoot/FF*1000,'k',xp,-chemPot_gmaShoot/FF*1000,'r');
xlabel('y in AgLi_y')
ylabel('Li equilibrium voltage (mV)')
legend('Overshoot','Thermo');

%% Save
% 1. Save the Li chemical potential for future usage, in unit J/mol
chemP_AgLi_db=[ypp,chemPot_thermo,chemPot_overShoot];
save('chemP_AgLi_db.mat','chemP_AgLi_db')
% % 2. Save as csv format
% eqmVdata=[ypp(1:9240),-chemPot_thermo(1:9240)/FF*1000];
% csvwrite('eqmV_AgLi.csv',eqmVdata)

%% Define material functions
% 1. The Gibbs free energy of AgLi vs. Li content (x)
function freEng = freEng_AgLi(x,AA,phaseID)
    freEng = AA(phaseID,1) + AA(phaseID,2)*x + AA(phaseID,3)*x.^2 + AA(phaseID,4)*x.^3 + ...
             AA(phaseID,5)*x.^4 + AA(phaseID,6)*x.^5 + AA(phaseID,7)*x.^6;      % In unit eV/atom
    freEng = freEng * 96487;                                                    % In unit J/mol. Unit covert: 1 eV/atom = 96487 J/mol
end

% 2. The Li chemical potential of AgLi_x vs. Li content (x), in unit J/mol
function chemPot = chemPot_AgLi(x,AA,phaseID)
    freEng  = freEng_AgLi(x,AA,phaseID);
    chemPot = freEng + 96487*(1-x) .* (AA(phaseID,2) + 2*AA(phaseID,3)*x + 3*AA(phaseID,4)*x.^2 + 4*AA(phaseID,5)*x.^3 + ...
                                       5*AA(phaseID,6)*x.^4 + 6*AA(phaseID,7)*x.^5);
end