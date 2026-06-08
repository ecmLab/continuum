% pulse.m
% ECRAM theoretical support -- single-electrode diffuse double layer under a
% PULSE TRAIN gate excitation. First-principles (COMSOL PNP) counterpart of the
% Python result fig1_accumulation.png: shows that a fast pulse train pre-loads
% the channel/electrolyte interface with Li+ (cumulative accumulation), whereas
% a single pulse fully relaxes (volatile EDL).
%
% Configuration (matches eis_ref.m single-electrode setup):
%   x=0  : blocking electrode (channel side) with Stern layer, driven by phiM(t)
%   x=L  : bulk reservoir (ground potential, fixed bulk concentration)
%
% Drive: phiM(t) = pulse train, amplitude V_amp, on-time t_on, off-time t_off.
%   Set N_PULSE and the on/off times to compare single vs cumulative pulsing.
%   Run once with t_off << relaxation (accumulation) and once with N_PULSE=1
%   (single pulse) to reproduce the volatile vs pre-loading contrast.
%
% Based on sinusoidal.m. Created for project_ecram: 2026-06-07

clc; clear; close all;

import com.comsol.model.*
import com.comsol.model.util.*

%% ---------------- Pulse-train configuration ----------------
V_dc    = 0.0;      % baseline electrode potential (V)
V_amp   = -0.20;    % pulse amplitude (V); NEGATIVE drives cation (Li+) to x=0
t_on    = 0.5;      % pulse ON duration (s)
t_off   = 0.2;      % pulse OFF duration (s)  -> set small for fast repetition
N_PULSE = 10;       % number of pulses (set 1 for the single-pulse control)
t_relax = 3.0;      % final relaxation window after the train (s)
t_rise  = 0.02*t_on;% smoothing of pulse edges (s)
pts_per_seg = 60;   % output points per on/off segment

period  = t_on + t_off;
t_total = N_PULSE*period + t_relax;

fprintf('=== ECRAM PULSE-TRAIN PNP ===\n');
fprintf('  V_amp=%.3f V, t_on=%.3f s, t_off=%.3f s, N=%d, t_total=%.2f s\n', ...
        V_amp, t_on, t_off, N_PULSE, t_total);

%% ---------------- Build pulse-train waveform table (MATLAB side) -----------
% smoothed square pulses via tanh edges; fed to COMSOL as an Interpolation fn
tt = linspace(0, t_total, 4000)';
VV = zeros(size(tt));
for k = 0:N_PULSE-1
    t0 = k*period;
    up   = 0.5*(1+tanh((tt - t0)/t_rise));
    down = 0.5*(1+tanh((t0 + t_on - tt)/t_rise));
    VV = VV + up.*down;
end
VV = V_dc + V_amp*min(VV,1);   % clamp overlap, apply amplitude

%% ---------------- COMSOL model ----------------
model = ModelUtil.create('Model');
model.modelPath(pwd);

model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 1);
model.component('comp1').mesh.create('mesh1');

% Electrostatics (phi) + Diluted Species (cA cation, cX anion), coupled (PNP)
model.component('comp1').physics.create('es', 'Electrostatics', 'geom1');
model.component('comp1').physics('es').field('electricpotential').field('phi');
model.component('comp1').physics.create('tds', 'DilutedSpecies', {'cA' 'cX'});
model.component('comp1').physics('tds').prop('TransportMechanism').set('Convection', '0');
model.component('comp1').physics('tds').prop('TransportMechanism').set('Migration', '1');
model.component('comp1').physics('tds').prop('ShapeProperty').set('order_concentration', '2');

model.component('comp1').multiphysics.create('pc1', 'PotentialCoupling', 1);
model.component('comp1').multiphysics('pc1').set('PotentialSource_physics', 'es');
model.component('comp1').multiphysics('pc1').set('PotentialDestination_physics', 'tds');
model.component('comp1').multiphysics('pc1').selection.all;
model.component('comp1').multiphysics.create('scdc1', 'SpaceChargeDensityCoupling', 1);
model.component('comp1').multiphysics('scdc1').set('SpaceChargeDensitySource_physics', 'tds');
model.component('comp1').multiphysics('scdc1').set('SpaceChargeDensityDestination_physics', 'es');
model.component('comp1').multiphysics('scdc1').selection.all;

model.study.create('std1');
model.study('std1').create('time', 'Transient');

%% ---------------- Parameters (solid Li-salt electrolyte) ----------------
model.param.set('T0', '25[degC]', 'Temperature');
model.param.set('V_therm', 'R_const*T0/F_const', 'Thermal voltage');
model.param.set('DA', '1e-13[m^2/s]', 'Diffusion coefficient, Li+ cation (slow, solid electrolyte)');
model.param.set('DX', 'DA*2', 'Diffusion coefficient, anion');
model.param.set('cA_bulk', '100[mol/m^3]', 'Bulk cation concentration');
model.param.set('cX_bulk', 'cA_bulk', 'Bulk anion concentration');
model.param.set('zA', '+1', 'Cation charge');
model.param.set('zX', '-1', 'Anion charge');
model.param.set('Istr_bulk', '0.5*((zA^2+zX^2)*cA_bulk)', 'Bulk ionic strength');
model.param.set('eps_PEO', '10.0', 'Relative permittivity of electrolyte');
model.param.set('eps_S', '10', 'Relative permittivity of Stern layer');
model.param.set('xD', 'sqrt(epsilon0_const*eps_PEO*V_therm/(2*F_const*Istr_bulk))', 'Debye length');
model.param.set('xS', '0.3[nm]', 'Stern layer thickness');
% Cell length chosen so tau_D = L^2/DA ~ a few s  (accumulation regime for sub-s pulses)
model.param.set('L_cell', '0.7[um]', 'Electrolyte thickness (tau_D = L^2/DA ~ 5 s)');
model.param.set('h_max', 'L_cell/40', 'Max mesh size (bulk)');
model.param.set('h_max_surf', 'xD/10', 'Max mesh size (electrode)');

% Pulse waveform as interpolation function pulseV(t)
model.func.create('pulseV', 'Interpolation');
model.func('pulseV').set('source', 'table');
model.func('pulseV').set('table', num2cell([tt VV]));   % [t(s) V(V)]
model.func('pulseV').set('interp', 'linear');
model.func('pulseV').set('extrap', 'const');
model.func('pulseV').set('argunit', 's');
model.func('pulseV').set('fununit', 'V');

model.param.set('phiM_R', '0[V]', 'Bulk reservoir potential (ground)');

%% ---------------- Stern-layer surface-charge variables ----------------
model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('phiM_L', 'pulseV(t)', 'Driven electrode potential (pulse train)');
model.component('comp1').variable('var1').set('deltaphi_L', 'phiM_L-phi', 'Electrode-OHP potential difference');
model.component('comp1').variable('var1').set('rho_surf_L', 'epsilon0_const*eps_S*deltaphi_L/xS', 'Stern surface charge density');
model.component('comp1').variable('var1').set('i_charging_L', 'd(rho_surf_L,t)', 'Charging current density (A/m^2)');
model.component('comp1').variable('var1').set('enrich_A', 'cA/cA_bulk', 'Cation enrichment factor');

%% ---------------- Geometry / BCs ----------------
model.component('comp1').geom('geom1').create('i1', 'Interval');
model.component('comp1').geom('geom1').feature('i1').setIndex('coord', 'L_cell', 1);
model.component('comp1').geom('geom1').run;

model.component('comp1').physics('es').feature('ccn1').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('es').feature('ccn1').set('epsilonr', {'eps_PEO' '0' '0' '0' 'eps_PEO' '0' '0' '0' 'eps_PEO'});

% x=0 : Stern surface charge (blocking electrode)
model.component('comp1').physics('es').create('sfcd1', 'SurfaceChargeDensity', 0);
model.component('comp1').physics('es').feature('sfcd1').selection.set([1]);
model.component('comp1').physics('es').feature('sfcd1').set('rhoqs', 'rho_surf_L');
% x=L : ground (bulk reservoir)
model.component('comp1').physics('es').create('pot1', 'ElectricPotential', 0);
model.component('comp1').physics('es').feature('pot1').selection.set([2]);
model.component('comp1').physics('es').feature('pot1').set('V0', 'phiM_R');

% species: x=0 blocking (no flux, default); x=L bulk concentration
model.component('comp1').physics('tds').feature('sp1').setIndex('z', 'zA', 0);
model.component('comp1').physics('tds').feature('sp1').setIndex('z', 'zX', 1);
model.component('comp1').physics('tds').feature('cdm1').set('D_cA', {'DA' '0' '0' '0' 'DA' '0' '0' '0' 'DA'});
model.component('comp1').physics('tds').feature('cdm1').set('D_cX', {'DX' '0' '0' '0' 'DX' '0' '0' '0' 'DX'});
model.component('comp1').physics('tds').feature('init1').setIndex('initc', 'cA_bulk', 0);
model.component('comp1').physics('tds').feature('init1').setIndex('initc', 'cX_bulk', 1);
model.component('comp1').physics('tds').create('conc1', 'Concentration', 0);
model.component('comp1').physics('tds').feature('conc1').selection.set([2]);
model.component('comp1').physics('tds').feature('conc1').setIndex('species', true, 0);
model.component('comp1').physics('tds').feature('conc1').setIndex('species', true, 1);
model.component('comp1').physics('tds').feature('conc1').setIndex('c0', 'cA_bulk', 0);
model.component('comp1').physics('tds').feature('conc1').setIndex('c0', 'cX_bulk', 1);

model.common('cminpt').set('modified', {'temperature' 'T0'});

%% ---------------- Mesh (refined at electrode) ----------------
model.component('comp1').mesh('mesh1').create('edg1', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg1').create('size_L', 'Size');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmax', 'h_max');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_L').selection.geom('geom1', 0);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_L').selection.set([1]);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_L').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_L').set('hmax', 'h_max_surf');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_L').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').run;

%% ---------------- Time list ----------------
% dense output across the whole train + relaxation
t_list = unique([linspace(0, t_total, N_PULSE*2*pts_per_seg + 200)]);
t_list_str = sprintf('%.10e ', t_list);
model.study('std1').feature('time').set('tlist', t_list_str);
% relative tolerance for the stiff EDL dynamics
model.study('std1').feature('time').set('rtol', '1e-4');

fprintf('Running transient pulse-train solve...\n'); tic;
model.study('std1').run();
fprintf('  done in %.1f s\n', toc);

%% ---------------- Extract interfacial accumulation ----------------
ti = mpheval(model, {'t','phiM_L','phi','cA','cX','enrich_A','rho_surf_L','i_charging_L'}, ...
             'edim', 0, 'selection', 1);
t_data   = ti.d1;   phiM_t = ti.d2;  phi0_t = ti.d3;
cA0_t    = ti.d4;   cX0_t  = ti.d5;  enr_t  = ti.d6;
rho_t    = ti.d7;   icharg_t = ti.d8;

%% ---------------- Output ----------------
out = 'rst/pulse';
if ~exist(out,'dir'); mkdir(out); end
mphsave(model, fullfile(out,'pulse.mph'));
save(fullfile(out,'pulse_results.mat'), 't_data','phiM_t','phi0_t','cA0_t', ...
     'cX0_t','enr_t','rho_t','icharg_t','V_amp','t_on','t_off','N_PULSE');
writetable(table(t_data, phiM_t, cA0_t, enr_t, rho_t, ...
    'VariableNames', {'t_s','phiM_V','cA_interface_molm3','enrichment','rho_surf_Cm2'}), ...
    fullfile(out,'pulse_interface.csv'));

figure;
subplot(2,1,1);
plot(t_data, phiM_t,'r-','LineWidth',1.5); ylabel('\phi_M (V)');
title(sprintf('Pulse train (N=%d, t_{on}=%.2fs, t_{off}=%.2fs)', N_PULSE, t_on, t_off)); grid on;
subplot(2,1,2);
plot(t_data, enr_t,'b-','LineWidth',1.5); hold on; yline(1,'k:');
xlabel('time (s)'); ylabel('interfacial c_A / c_{bulk}');
title('Interfacial Li^+ enrichment (accumulation under fast pulsing)'); grid on;
saveas(gcf, fullfile(out,'pulse_accumulation.png'));

fprintf('\nInterfacial enrichment: peak=%.2f x bulk, final=%.2f x bulk\n', ...
        max(enr_t), enr_t(end));
fprintf('Results in ./%s/  (compare with python fig1_accumulation.png)\n', out);
fprintf('TIP: rerun with N_PULSE=1 for the single-pulse (volatile) control.\n');
