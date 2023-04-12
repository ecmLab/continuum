% clc;clear;
% set(0, 'DefaultLineLineWidth', 2);
% set(0,'defaultAxesFontSize',16);
AA    = [0.1825, -61.39, 59.00, 0.0, 0.0, 0.0; 1916.9, -14062,  40108, -56284, 38879, -10544;...
         29.53,  -164.4, 141.1, 0.0, 0.0, 0.0; -51.169, 50.806, 0.0,   0.0,    0.0,   0.0];
tang_xp = [0.16798118, 0.4894171, 0.568393, 0.7102594];   % The tangent points of the Gibbs free energy between phases
tang_yp = tang_xp ./ (1-tang_xp);                        % The tangent points of the equilibrium potential between phases
Npp     = 10000;

% The free energy of three phases, in unit J/mol
xp          = linspace(0.0,0.99,Npp)';
freEng_fcc  = freEng_AgLi(xp,AA,1);
freEng_fcc(xp>0.75) = NaN;
freEng_beta = freEng_AgLi(xp,AA,2);
freEng_beta(xp<0.35) = NaN;
freEng_gma  = freEng_AgLi(xp,AA,3);
freEng_gma(xp<0.4) = NaN;
freEng_bcc  = freEng_AgLi(xp,AA,4);
freEng_bcc(xp<0.75) = NaN;

% The Li chemical potential of each phase, in unit J/mol
chemPot_fcc = chemPot_AgLi(xp,AA,1);
chemPot_beta= chemPot_AgLi(xp,AA,2);
chemPot_gma = chemPot_AgLi(xp,AA,3);
chemPot_bcc = chemPot_AgLi(xp,AA,4);
chemPot_gma(chemPot_gma>chemPot_bcc(end)) = chemPot_bcc(end);

% The Gibbs free energy and Li chemical potential of the AgLiy with overshoot
freEng_overShoot  = zeros(Npp,1);
chemPot_overShoot = zeros(Npp,1);
for i = 1 : Npp
    [freEng_overShoot(i),idx] = min([freEng_fcc(i),freEng_beta(i),freEng_gma(i),freEng_bcc(i)]);
    eqm = [chemPot_fcc(i),chemPot_beta(i),chemPot_gma(i),chemPot_bcc(i)];
    chemPot_overShoot(i) = eqm(idx);
end

% The Gibbs free energy and Li chemical potential of the AgLiy at thermodynamic equilibrium
freEng_thermo  = freEng_overShoot;
chemPot_thermo = chemPot_overShoot;
% Find the co-tangent line between FCC and beta phase
idx_fcc_beta                  = find(xp>tang_xp(1) & xp<tang_xp(2));
chemPot_thermo(idx_fcc_beta)  = chemPot_thermo(idx_fcc_beta(1)-1);
dedx                          = (freEng_thermo(idx_fcc_beta(end))-freEng_thermo(idx_fcc_beta(1)))/(xp(idx_fcc_beta(end))-xp(idx_fcc_beta(1)));
freEng_thermo(idx_fcc_beta)   = dedx * (xp(idx_fcc_beta)-xp(idx_fcc_beta(1))) + freEng_thermo(idx_fcc_beta(1));
% Find the co-tangent line between beta and gamma phase
idx_bcc_gma                   = find(xp>tang_xp(3) & xp<tang_xp(4));
chemPot_thermo(idx_bcc_gma)   = chemPot_thermo(idx_bcc_gma(1)-1);
dedx                          = (freEng_thermo(idx_bcc_gma(end))-freEng_thermo(idx_bcc_gma(1)))/(xp(idx_bcc_gma(end))-xp(idx_bcc_gma(1)));
freEng_thermo(idx_bcc_gma)    = dedx * (xp(idx_bcc_gma)-xp(idx_bcc_gma(1))) + freEng_thermo(idx_bcc_gma(1));
% % Find the co-tangent line between gamma1 and gamma2 phase
% idx_gma1_gma2                = find(xp>tang_xp(5) & xp<tang_xp(6));
% chemPot_thermo(idx_gma1_gma2)= chemPot_thermo(idx_gma1_gma2(2)-1);
% dedx                         = (freEng_thermo(idx_gma1_gma2(end))-freEng_thermo(idx_gma1_gma2(1)))/(xp(idx_gma1_gma2(end))-xp(idx_gma1_gma2(1)));
% freEng_thermo(idx_gma1_gma2) = dedx * (xp(idx_gma1_gma2)-xp(idx_gma1_gma2(1))) + freEng_thermo(idx_gma1_gma2(1));
% % Find the co-tangent line between gamma2 and bcc Li phase
% idx_gma2_bcc                 = find(xp>tang_xp(7) & xp<tang_xp(8));
% chemPot_thermo(idx_gma2_bcc) = chemPot_thermo(idx_gma2_bcc(1)-1);
% dedx                         = (freEng_thermo(idx_gma2_bcc(end))-freEng_thermo(idx_gma2_bcc(1)))/(xp(idx_gma2_bcc(end))-xp(idx_gma2_bcc(1)));
% freEng_thermo(idx_gma2_bcc)  = dedx * (xp(idx_gma2_bcc)-xp(idx_gma2_bcc(1))) + freEng_thermo(idx_gma2_bcc(1));

% The Li chemical potential of the AgLiy when voltage overshoot to zero from BCC phase;
chemPot_betaShoot = chemPot_thermo;
for i = 1 : Npp
    if xp(i) > tang_xp(4)
        chemPot_betaShoot(i) = chemPot_beta(i);
    end
end
% idtp = find(chemPot_betaShoot>0);
% chemPot_betaShoot(idtp) = 0;

% % The Li chemical potential of the AgLiy when voltage overshoot to zero from gama phase;
% chemPot_gmaShoot = chemPot_thermo;
% for i = 1 : Npp
%     if xp(i) > tang_xp(7)
%         chemPot_gmaShoot(i) = chemPot_gma(i);
%     end
% end
% idtp = find(chemPot_gmaShoot>0);
% chemPot_gmaShoot(idtp) = 0;

%% Plot
% 1. Plot the free energy of each phase
ifg = ifg + 1;
figure(ifg)
plot(xp,freEng_fcc,'-.',xp,freEng_beta,'-.',xp,freEng_gma,'--',xp,freEng_bcc,'--');
xlabel('x in Ag_{1-x}Li_x')
ylabel('Gibbs free energy (J/mol)')
legend('FCC','Beta','GAMMA','BCC');
% 2. Plot the overall free energy and equilibrium voltage vs. x in Ag_{1-x}Li_x
ifg = ifg + 1;
figure(ifg)
plot(xp,freEng_overShoot,xp,freEng_thermo,'--');
xlabel('x in Ag_{1-x}Li_x')
ylabel('Gibbs free energy (J/mol)')
legend('Overshoot','Thermo');
ifg = ifg + 1;
figure(ifg)
plot(xp,chemPot_overShoot,'r',xp,chemPot_betaShoot,'b',xp,chemPot_thermo,'--k');
xlabel('x in Ag_{1-x}Li_x')
ylabel('Li chemical potential (J/mol)')
legend('Overshoot','betaShoot','Thermo');

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
plot(ypp,-chemPot_overShoot/FF*1000,'r',ypp,-chemPot_betaShoot/FF*1000,'b',ypp,-chemPot_thermo/FF*1000,'--k');
xlabel('y in AgLi_y')
ylabel('Li equilibrium voltage (mV)')
legend('Overshoot','betaShoot','Thermo');

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
    freEng = AA(phaseID,1) + AA(phaseID,2)*x + AA(phaseID,3)*x.^2 + AA(phaseID,4)*x.^3 + AA(phaseID,5)*x.^4 + AA(phaseID,6)*x.^5; % In unit kJ/mol
    freEng = freEng * 1000;                                  % Convert to unit J/mol
end

% 2. The Li chemical potential of AgLi_x vs. Li content (x), in unit J/mol
function chemPot = chemPot_AgLi(x,AA,phaseID)
    freEng  = freEng_AgLi(x,AA,phaseID);
    chemPot = freEng + 1000*(1-x) .* (AA(phaseID,2) + 2*AA(phaseID,3)*x + 3*AA(phaseID,4)*x.^2 + 4*AA(phaseID,5)*x.^3 + 5*AA(phaseID,6)*x.^4);
end