% clc;clear;
% set(0, 'DefaultLineLineWidth', 2);
% set(0,'defaultAxesFontSize',16);

AA      = [-0.00118696639, -1.19634947,  4.82894652, -18.1038119,  35.9675069,  -32.2379316,  10.7413637;...
           -0.03263715,    -0.98628971,  3.36536168, -11.34600958, 21.8986608,  -19.24766004, 6.34098313;...
            516.1876733,    -3922.736263, 12356.5479, -20663.68228, 19349.86321, -9620.179228, 1984.004148];
tang_xp = [0.2965,0.3639,0.5552,0.6411,0.8219,0.8846];   % The tangent points of the Gibbs free energy between phases
tang_yp = tang_xp ./ (1-tang_xp);                        % The tangent points of the equilibrium potential between phases
Npp     = 10000;

% The free energy of three phases, in unit J/mol
xp          = linspace(0.0,yp_end/(1+yp_end),Npp)';
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
    [freEng_overShoot(i),idx] = min([freEng_bcc(i),freEng_fcc(i),freEng_gma(i)]);
    eqm = [chemPot_bcc(i),chemPot_fcc(i),chemPot_gma(i)];
    chemPot_overShoot(i) = eqm(idx);
end

% The Gibbs free energy and Li chemical potential of the AgLiy at thermodynamic equilibrium
freEng_thermo  = freEng_overShoot;
chemPot_thermo = chemPot_overShoot;
% Find the co-tangent line between bcc and fcc phase
idx_bcc_fcc    = find(xp>tang_xp(1) & xp<tang_xp(2));
chemPot_thermo(idx_bcc_fcc) = chemPot_thermo(idx_bcc_fcc(1)-1);
dedx  = (freEng_thermo(idx_bcc_fcc(end))-freEng_thermo(idx_bcc_fcc(1)))/(xp(idx_bcc_fcc(end))-xp(idx_bcc_fcc(1)));
freEng_thermo(idx_bcc_fcc) = dedx * (xp(idx_bcc_fcc)-xp(idx_bcc_fcc(1))) + freEng_thermo(idx_bcc_fcc(1));
% Find the co-tangent line between bcc and gamma left phase
idx_bcc_gma_left = find(xp>tang_xp(3) & xp<tang_xp(4));
chemPot_thermo(idx_bcc_gma_left) = chemPot_thermo(idx_bcc_gma_left(1)-1);
dedx  = (freEng_thermo(idx_bcc_gma_left(end))-freEng_thermo(idx_bcc_gma_left(1)))/(xp(idx_bcc_gma_left(end))-xp(idx_bcc_gma_left(1)));
freEng_thermo(idx_bcc_gma_left) = dedx * (xp(idx_bcc_gma_left)-xp(idx_bcc_gma_left(1))) + freEng_thermo(idx_bcc_gma_left(1));
% Find the co-tangent line between bcc and gamma right phase
idx_bcc_gma_rght = find(xp>tang_xp(5) & xp<tang_xp(6));
chemPot_thermo(idx_bcc_gma_rght) = chemPot_thermo(idx_bcc_gma_rght(1)-1);
dedx  = (freEng_thermo(idx_bcc_gma_rght(end))-freEng_thermo(idx_bcc_gma_rght(1)))/(xp(idx_bcc_gma_rght(end))-xp(idx_bcc_gma_rght(1)));
freEng_thermo(idx_bcc_gma_rght) = dedx * (xp(idx_bcc_gma_rght)-xp(idx_bcc_gma_rght(1))) + freEng_thermo(idx_bcc_gma_rght(1));

%% Plot
% 1. Plot the free energy of each phase
ifg = ifg + 1;
figure(ifg)
plot(xp,freEng_bcc,'-.',xp,freEng_fcc,'-.',xp(Npp*0.62:end),freEng_gma(Npp*0.62:end),'--');
xlabel('x in Ag_{1-x}Li_x')
ylabel('Gibbs free energy (J/mol)')
legend('BCC','FCC','GAMMA');
% 2. Plot the overall free energy and equilibrium voltage vs. x in Ag_{1-x}Li_x
ifg = ifg + 1;
figure(ifg)
% subplot(2,1,1);
plot(xp,freEng_overShoot,xp,freEng_thermo,'--');
xlabel('x in Ag_{1-x}Li_x')
ylabel('Gibbs free energy (J/mol)')
legend('Overshoot','Thermo');
ifg = ifg + 1;
figure(ifg)
% subplot(2,1,2);
plot(xp,chemPot_overShoot,xp,chemPot_thermo,'--');
xlabel('x in Ag_{1-x}Li_x')
ylabel('Li chemical potential (J/mol)')
legend('Overshoot','Thermo');

% 3. Plot the overall free energy and equilibrium voltage vs. y in AgLi_y
ifg = ifg + 1;
figure(ifg)
ypp = xp ./ (1-xp);
% subplot(2,1,1);
plot(ypp,freEng_overShoot./(1-xp),ypp,freEng_thermo./(1-xp),'--');
xlabel('y in AgLi_y')
ylabel('Gibbs free energy (J/mol)')
legend('Overshoot','Thermo');
ifg = ifg + 1;
figure(ifg)
ypp = xp ./ (1-xp);
% subplot(2,1,2);
% plot(yp,chemPot_bcc,'--',yp,chemPot_fcc,'--',yp(Np*0.65:end),chemPot_gma(Np*0.65:end),'--',yp,chemPot_overShoot,yp,chemPot_thermo);
plot(ypp,chemPot_overShoot,ypp,chemPot_thermo,'--');
xlabel('y in AgLi_y')
ylabel('Li chemical potential (J/mol)')
legend('Overshoot','Thermo');

%% Save
% 1. Save the Li chemical potential for future usage, in unit J/mol
chemP_AgLi_db=[ypp/yp_end,chemPot_thermo,chemPot_overShoot];
save('chemP_AgLi_db.mat','chemP_AgLi_db')
% % 2. Save as csv format
% eqmVdata=[ypp,chemPot_thermo*1000,chemPot_overShoot*1000];
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