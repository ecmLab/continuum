% twoInterface_v1.m
% Two-interface diffuse double layer model - Two-Electrode Device
% Left electrode at V with Stern layer, Right electrode GROUNDED with Stern layer
% Created: 2025-10-09
% Updated: 2025-10-20 - Parameters aligned with Mei et al. 2018 JPCC paper
%
% Reference: Mei, B.-A., Munteshari, O., Lau, J., Dunn, B., and Pilon, L.
%            "Physical Interpretations of Nyquist Plots for EDLC Electrodes and Devices"
%            J. Phys. Chem. C 2018, 122, 194-206
%            DOI: 10.1021/acs.jpcc.7b10582

clc; clear;
close all;

%% Setup COMSOL
import com.comsol.model.*
import com.comsol.model.util.*

fprintf('Building TWO-ELECTRODE DEVICE model...\n');
fprintf('Configuration: Left electrode at V (with Stern layer), Right electrode grounded (with Stern layer)\n');
fprintf('Parameters aligned with Mei et al. 2018 JPCC paper\n\n');

% Create new model
model = ModelUtil.create('Model');
model.modelPath(pwd);

% Create component and geometry (1D)
model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 1);
model.component('comp1').mesh.create('mesh1');

% Create physics: Electrostatics
model.component('comp1').physics.create('es', 'Electrostatics', 'geom1');
model.component('comp1').physics('es').field('electricpotential').field('phi');

% Create physics: Diluted Species (cation and anion)
model.component('comp1').physics.create('tds', 'DilutedSpecies', {'cA' 'cX'});
model.component('comp1').physics('tds').prop('TransportMechanism').set('Convection', '0');
model.component('comp1').physics('tds').prop('TransportMechanism').set('Migration', '1');
model.component('comp1').physics('tds').prop('ShapeProperty').set('order_concentration', '2');

% Create multiphysics couplings
model.component('comp1').multiphysics.create('pc1', 'PotentialCoupling', 1);
model.component('comp1').multiphysics('pc1').set('PotentialSource_physics', 'es');
model.component('comp1').multiphysics('pc1').set('PotentialDestination_physics', 'tds');
model.component('comp1').multiphysics('pc1').selection.all;

model.component('comp1').multiphysics.create('scdc1', 'SpaceChargeDensityCoupling', 1);
model.component('comp1').multiphysics('scdc1').set('SpaceChargeDensitySource_physics', 'tds');
model.component('comp1').multiphysics('scdc1').set('SpaceChargeDensityDestination_physics', 'es');
model.component('comp1').multiphysics('scdc1').selection.all;

% Create study with both stationary and time-dependent steps
model.study.create('std1');
model.study('std1').create('stat', 'Stationary');
model.study('std1').create('time', 'Transient');

% Set parameters (aligned with Mei et al. 2018 JPCC paper and oneInterface_v1.m)
% Reference: Mei, B.-A. et al. J. Phys. Chem. C 2018, 122, 194-206
% DOI: 10.1021/acs.jpcc.7b10582

% Temperature and fundamental constants
model.param.set('T0', '298[K]', 'Temperature (298 K from paper)');
model.param.set('V_therm', 'R_const*T0/F_const', 'Thermal voltage');

% Ion transport properties (Table 1 from paper)
model.param.set('DA', '2e-13[m^2/s]', 'Diffusion coefficient, cation (typical value from paper)');
model.param.set('DX', 'DA', 'Diffusion coefficient, anion (symmetric electrolyte)');
model.param.set('cA_bulk', '1[mol/m^3]', 'Bulk cation concentration (0.001 mol/L from paper Table 2)');
model.param.set('cX_bulk', 'cA_bulk', 'Bulk anion concentration');
model.param.set('zA', '+1', 'Cation charge (z=1 from paper)');
model.param.set('zX', '-1', 'Anion charge');
model.param.set('Istr_bulk', '0.5*((zA^2+zX^2)*cA_bulk)', 'Bulk ionic strength');

% Electrolyte properties (Table 1 from paper)
model.param.set('eps_r', '64.4', 'Relative permittivity of electrolyte (propylene carbonate from paper)');
model.param.set('a_ion', '0.66[nm]', 'Effective ion diameter (from paper Table 2)');
model.param.set('xS', 'a_ion/2', 'Stern layer thickness (H = a/2 from paper)');
model.param.set('xD', 'sqrt(epsilon0_const*eps_r*V_therm/(2*F_const*Istr_bulk))', 'Debye length');

% Domain dimensions (Table 1 from paper)
% For two-electrode device: total length = 2×L_electrode (Case 25 in Table 2: 2×1600nm = 3200nm)
% Using smaller domain for faster computation: 2×160nm = 320nm
model.param.set('L_cell', '320[nm]', 'Total cell length (2 electrodes separated by electrolyte)');
model.param.set('Le', '10[nm]', 'Electrode thickness (from paper Table 1)');

% Mesh parameters
model.param.set('h_max', 'L_cell/80', 'Maximum mesh element size (finer for two electrodes)');
model.param.set('h_max_surf', 'xD/100', 'Maximum mesh element size (electrodes)');

% Capacitance (for reference)
model.param.set('Cd_GCS', 'epsilon0_const/(xD/eps_r+xS/eps_r)', 'Capacitance per unit area (GCS theory)');

% Applied potential - TWO-ELECTRODE DEVICE CONFIGURATION
% Left electrode: Applied voltage, Right electrode: Grounded
model.param.set('phiM_L', '0.02[V]', 'Left electrode potential (applied voltage)');
model.param.set('phiM_R', '0[V]', 'Right electrode potential (grounded)');

% Create variables for BOTH electrodes (both have Stern layers!)
model.component('comp1').variable.create('var1');
% Left electrode variables
model.component('comp1').variable('var1').set('deltaphi_L', 'phiM_L-phi', 'Left electrode-OHP potential difference');
model.component('comp1').variable('var1').set('rho_surf_L', 'epsilon0_const*eps_r*deltaphi_L/xS', 'Left surface charge density');
model.component('comp1').variable('var1').set('i_charging_L', 'd(rho_surf_L,t)', 'Left charging current density (A/m^2)');
% Right electrode variables (grounded, has Stern layer)
model.component('comp1').variable('var1').set('deltaphi_R', 'phiM_R-phi', 'Right electrode-OHP potential difference');
model.component('comp1').variable('var1').set('rho_surf_R', 'epsilon0_const*eps_r*deltaphi_R/xS', 'Right surface charge density');
model.component('comp1').variable('var1').set('i_charging_R', 'd(rho_surf_R,t)', 'Right charging current density (A/m^2)');
% Flux and current variables
% Use tds.tflux (total flux) directly - it's the correct volume variable for bulk evaluation
model.component('comp1').variable('var1').set('i_cation', 'F_const*zA*tds.tflux_cAx', 'Cation current density (A/m^2)');
model.component('comp1').variable('var1').set('i_anion', 'F_const*zX*tds.tflux_cXx', 'Anion current density (A/m^2)');
model.component('comp1').variable('var1').set('i_total', 'i_cation+i_anion', 'Total ionic current density (A/m^2)');
% Charge conservation check
model.component('comp1').variable('var1').set('Q_total', 'rho_surf_L+rho_surf_R', 'Total charge (should be ~0 for capacitor)');

% Create geometry (1D interval from 0 to L_cell)
model.component('comp1').geom('geom1').create('i1', 'Interval');
model.component('comp1').geom('geom1').feature('i1').setIndex('coord', 'L_cell', 1);
model.component('comp1').geom('geom1').run;

% Set physics properties - Electrostatics
model.component('comp1').physics('es').feature('ccn1').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('es').feature('ccn1').set('epsilonr', {'eps_r' '0' '0' '0' 'eps_r' '0' '0' '0' 'eps_r'});

% Boundary conditions - Electrostatics
% Left electrode (boundary 1, x=0): Surface charge from Stern layer
model.component('comp1').physics('es').create('sfcd1', 'SurfaceChargeDensity', 0);
model.component('comp1').physics('es').feature('sfcd1').selection.set([1]);
model.component('comp1').physics('es').feature('sfcd1').set('rhoqs', 'rho_surf_L');

% Right electrode (boundary 2, x=L): Surface charge from Stern layer
model.component('comp1').physics('es').create('sfcd2', 'SurfaceChargeDensity', 0);
model.component('comp1').physics('es').feature('sfcd2').selection.set([2]);
model.component('comp1').physics('es').feature('sfcd2').set('rhoqs', 'rho_surf_R');

% NOTE: NO explicit ground condition - the grounded metal is enforced through phiM_R = 0
% Both electrodes have Stern layers, making this different from oneInterface

% Set physics properties - Transport of Diluted Species
model.component('comp1').physics('tds').feature('sp1').setIndex('z', 'zA', 0);
model.component('comp1').physics('tds').feature('sp1').setIndex('z', 'zX', 1);
model.component('comp1').physics('tds').feature('cdm1').set('D_cA', {'DA' '0' '0' '0' 'DA' '0' '0' '0' 'DA'});
model.component('comp1').physics('tds').feature('cdm1').set('D_cX', {'DX' '0' '0' '0' 'DX' '0' '0' '0' 'DX'});
model.component('comp1').physics('tds').feature('init1').setIndex('initc', 'cA_bulk', 0);
model.component('comp1').physics('tds').feature('init1').setIndex('initc', 'cX_bulk', 1);

% Transport BCs: Zero flux at both electrodes (blocking electrodes - default)

% Set common properties
model.common('cminpt').set('modified', {'temperature' 'T0'});

% Create mesh with refinement at BOTH electrodes
model.component('comp1').mesh('mesh1').create('edg1', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg1').create('size1', 'Size');  % General mesh
model.component('comp1').mesh('mesh1').feature('edg1').create('size_L', 'Size'); % Left electrode
model.component('comp1').mesh('mesh1').feature('edg1').create('size_R', 'Size'); % Right electrode

% General mesh
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmax', 'h_max');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmaxactive', true);

% Left electrode refinement (boundary 1)
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_L').selection.geom('geom1', 0);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_L').selection.set([1]);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_L').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_L').set('hmax', 'h_max_surf');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_L').set('hmaxactive', true);

% Right electrode refinement (boundary 2)
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_R').selection.geom('geom1', 0);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_R').selection.set([2]);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_R').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_R').set('hmax', 'h_max_surf');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_R').set('hmaxactive', true);

model.component('comp1').mesh('mesh1').run;

% Configure study settings
model.study('std1').feature('stat').active(false);

% Time-dependent study settings
model.study('std1').feature('time').set('tunit', 'ms');
% Use finer steps under 1 ms while keeping the original span to 1000 ms (double resolution <0.01 ms)
model.study('std1').feature('time').set('tlist', '10^{range(log10(0.001),1/12,log10(0.01))} 10^{range(log10(0.01)+1/6,1/6,log10(1))} 10^{range(log10(1),1/3,log10(1000))}');

fprintf('Two-interface (asymmetric) model built successfully!\n');
fprintf('  Model tag: %s\n', char(model.tag()));
fprintf('  Configuration: Asymmetric - V + grounded (both with Stern)\n');
fprintf('  Left electrode: +%.1f mV (with Stern layer)\n', 10.0);
fprintf('  Right electrode: 0 mV (grounded metal with Stern layer)\n');
fprintf('  Domain length: 200 x Debye length\n\n');

%% Run the time-dependent study
fprintf('Running time-dependent simulation...\n');
fprintf('Time range: 0.001 to 1000 ms (logarithmic)\n');
tic;
model.study('std1').run();
elapsed_time = toc;
fprintf('Simulation completed in %.2f seconds\n', elapsed_time);

%% Extract results along the domain (final time step)
fprintf('\nExtracting spatial profiles at final time...\n');

eval_result = mpheval(model, {'x', 'phi', 'cA', 'cX', 'es.normE', 'es.rhoq'}, 'edim', 1, 'solnum', 'end');

x_data = eval_result.d1;
phi_data = eval_result.d2;
cA_data = eval_result.d3;
cX_data = eval_result.d4;
Es_data = eval_result.d5;
rho_data = eval_result.d6;

fprintf('  Number of spatial points: %d\n', length(x_data));

%% Extract time series data at BOTH electrodes
fprintf('\nExtracting time evolution at left electrode (x=0)...\n');
time_L = mpheval(model, {'t', 'phi', 'cA', 'cX', 'i_charging_L', 'i_cation', 'i_anion', 'i_total', 'rho_surf_L'}, ...
    'edim', 0, 'selection', 1);
t_data = time_L.d1;
phi_L_t = time_L.d2;
cA_L_t = time_L.d3;
cX_L_t = time_L.d4;
i_charging_L_t = time_L.d5;
i_cation_L_t = time_L.d6;
i_anion_L_t = time_L.d7;
i_total_L_t = time_L.d8;
rho_L_t = time_L.d9;

fprintf('  Number of time points: %d\n', length(t_data));

fprintf('\nExtracting time evolution at right electrode (x=L)...\n');
time_R = mpheval(model, {'phi', 'cA', 'cX', 'i_charging_R', 'i_cation', 'i_anion', 'i_total', 'rho_surf_R'}, ...
    'edim', 0, 'selection', 2);
phi_R_t = time_R.d1;
cA_R_t = time_R.d2;
cX_R_t = time_R.d3;
i_charging_R_t = time_R.d4;
i_cation_R_t = time_R.d5;
i_anion_R_t = time_R.d6;
i_total_R_t = time_R.d7;
rho_R_t = time_R.d8;

fprintf('\nExtracting time evolution at cell midpoint (x=L/2)...\n');
% Get the actual midpoint coordinate value
x_mid = mphglobal(model, 'L_cell/2');
x_mid = x_mid(1); % Extract scalar value
fprintf('  Midpoint coordinate: %.4e m\n', x_mid);

% Use mpheval with edim=1 and recover='ppr' to get interior point values
% Then find the point closest to midpoint
num_times = length(t_data);
i_total_mid_t = zeros(num_times, 1);
i_cation_mid_t = zeros(num_times, 1);
i_anion_mid_t = zeros(num_times, 1);

% Get x coordinates from first solution to find midpoint index
eval_x = mpheval(model, 'x', 'edim', 1, 'solnum', 1);
x_coords = eval_x.d1;
[~, mid_idx] = min(abs(x_coords - x_mid)); % Find closest point to midpoint

fprintf('  Using mesh point at x = %.4e m (index %d)\n', x_coords(mid_idx), mid_idx);

% Extract current data at midpoint for all time steps
for idx = 1:num_times
    eval_i = mpheval(model, {'i_total', 'i_cation', 'i_anion'}, 'edim', 1, 'solnum', idx);
    i_total_mid_t(idx) = eval_i.d1(mid_idx);
    i_cation_mid_t(idx) = eval_i.d2(mid_idx);
    i_anion_mid_t(idx) = eval_i.d3(mid_idx);
end

% Calculate total charge for conservation check
Q_total_t = rho_L_t + rho_R_t;

%% Create output directory
output_dir = 'rst/twoInterface';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('\nCreated output directory: %s\n', output_dir);
end

%% Save COMSOL model file
model_file = fullfile(output_dir, 'twoInterface.mph');
mphsave(model, model_file);
fprintf('COMSOL model saved to: %s\n', model_file);

%% Save results to file
results_file = fullfile(output_dir, 'twoInterface_results.mat');
save(results_file, 'x_data', 'phi_data', 'cA_data', 'cX_data', 'Es_data', 'rho_data', ...
    't_data', 'phi_L_t', 'cA_L_t', 'cX_L_t', 'phi_R_t', 'cA_R_t', 'cX_R_t', ...
    'i_charging_L_t', 'i_cation_L_t', 'i_anion_L_t', 'i_total_L_t', 'rho_L_t', ...
    'i_charging_R_t', 'i_cation_R_t', 'i_anion_R_t', 'i_total_R_t', 'rho_R_t', ...
    'x_mid', 'i_cation_mid_t', 'i_anion_mid_t', 'i_total_mid_t', ...
    'Q_total_t', 'elapsed_time');
fprintf('Results saved to: %s\n', results_file);

%% Save spatial data as CSV
csv_spatial = fullfile(output_dir, 'twoInterface_spatial.csv');
spatial_table = table(x_data*1e9, phi_data*1000, cA_data, cX_data, Es_data/1e6, rho_data, ...
    'VariableNames', {'x_nm', 'phi_mV', 'cA_molm3', 'cX_molm3', 'E_MVm', 'rho_Cm3'});
writetable(spatial_table, csv_spatial);
fprintf('Spatial data saved to: %s\n', csv_spatial);

%% Save time evolution data
csv_time_L = fullfile(output_dir, 'twoInterface_time_left.csv');
time_L_table = table(t_data*1000, phi_L_t*1000, cA_L_t, cX_L_t, ...
    i_charging_L_t, i_cation_L_t, i_anion_L_t, i_total_L_t, rho_L_t, ...
    'VariableNames', {'t_ms', 'phi_mV', 'cA_molm3', 'cX_molm3', ...
    'i_charging_Am2', 'i_cation_Am2', 'i_anion_Am2', 'i_total_Am2', 'rho_surf_Cm2'});
writetable(time_L_table, csv_time_L);
fprintf('Left electrode time data saved to: %s\n', csv_time_L);

csv_time_R = fullfile(output_dir, 'twoInterface_time_right.csv');
time_R_table = table(t_data*1000, phi_R_t*1000, cA_R_t, cX_R_t, ...
    i_charging_R_t, i_cation_R_t, i_anion_R_t, i_total_R_t, rho_R_t, ...
    'VariableNames', {'t_ms', 'phi_mV', 'cA_molm3', 'cX_molm3', ...
    'i_charging_Am2', 'i_cation_Am2', 'i_anion_Am2', 'i_total_Am2', 'rho_surf_Cm2'});
writetable(time_R_table, csv_time_R);
fprintf('Right electrode time data saved to: %s\n', csv_time_R);

csv_time_mid = fullfile(output_dir, 'twoInterface_time_mid.csv');
time_mid_table = table(t_data*1000, i_total_mid_t, i_cation_mid_t, i_anion_mid_t, ...
    'VariableNames', {'t_ms', 'i_total_Am2', 'i_cation_Am2', 'i_anion_Am2'});
writetable(time_mid_table, csv_time_mid);
fprintf('Midpoint time data saved to: %s\n', csv_time_mid);

%% Create plots
fprintf('\nGenerating plots...\n');

% Figure 1: Spatial profiles (final state)
figure(1);
subplot(2,2,1);
plot(x_data*1e9, phi_data*1000, 'b-', 'LineWidth', 2);
xlabel('Position (nm)');
ylabel('Potential (mV)');
title('Electric Potential (Final State)');
grid on;

subplot(2,2,2);
plot(x_data*1e9, cA_data, 'r-', 'LineWidth', 2); hold on;
plot(x_data*1e9, cX_data, 'b-', 'LineWidth', 2);
xlabel('Position (nm)');
ylabel('Concentration (mol/m^3)');
legend('Cation A^+', 'Anion X^-');
title('Ion Concentrations (Final State)');
grid on;

subplot(2,2,3);
plot(x_data*1e9, Es_data/1e6, 'g-', 'LineWidth', 2);
xlabel('Position (nm)');
ylabel('E-field (MV/m)');
title('Electric Field (Final State)');
grid on;

subplot(2,2,4);
plot(x_data*1e9, rho_data, 'm-', 'LineWidth', 2);
xlabel('Position (nm)');
ylabel('Charge density (C/m^3)');
title('Space Charge Density (Final State)');
grid on;

saveas(gcf, fullfile(output_dir, 'twoInterface_spatial.png'));

% Figure 2: Time evolution comparison - Left vs Right
figure(2);
subplot(3,1,1);
semilogx(t_data*1000, phi_L_t*1000, 'r-', 'LineWidth', 2); hold on;
semilogx(t_data*1000, phi_R_t*1000, 'b-', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Potential (mV)');
title('Electrode Potentials vs Time');
legend('Left (V)', 'Right (ground)', 'Location', 'best');
grid on;

subplot(3,1,2);
semilogx(t_data*1000, cA_L_t, 'r-', 'LineWidth', 2); hold on;
semilogx(t_data*1000, cA_R_t, 'r--', 'LineWidth', 2);
semilogx(t_data*1000, cX_L_t, 'b-', 'LineWidth', 2);
semilogx(t_data*1000, cX_R_t, 'b--', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Concentration (mol/m^3)');
title('Ion Concentrations at Electrodes');
legend('A^+ Left', 'A^+ Right', 'X^- Left', 'X^- Right', 'Location', 'best');
grid on;

subplot(3,1,3);
semilogx(t_data*1000, i_total_L_t, 'r-', 'LineWidth', 2); hold on;
semilogx(t_data*1000, i_total_R_t, 'b-', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Current (A/m^2)');
title('Ionic Current at Electrodes');
legend('Left', 'Right', 'Location', 'best');
grid on;

saveas(gcf, fullfile(output_dir, 'twoInterface_time_comparison.png'));

% Figure 3: Charge conservation check
figure(3);
semilogx(t_data*1000, rho_L_t, 'r-', 'LineWidth', 2); hold on;
semilogx(t_data*1000, rho_R_t, 'b-', 'LineWidth', 2);
semilogx(t_data*1000, Q_total_t, 'k--', 'LineWidth', 2);
semilogx(t_data*1000, zeros(size(t_data)), 'k:', 'LineWidth', 1);
xlabel('Time (ms)');
ylabel('Surface Charge Density (C/m^2)');
title('Charge Conservation: Q_{left} + Q_{right} = 0');
legend('Q_{left}', 'Q_{right}', 'Q_{total}', 'Zero', 'Location', 'best');
grid on;
saveas(gcf, fullfile(output_dir, 'twoInterface_charge_conservation.png'));

% Figure 4: Current components at cell midpoint
figure(4);
semilogx(t_data*1000, i_total_mid_t, 'r-', 'LineWidth', 2); hold on;
semilogx(t_data*1000, i_cation_mid_t, 'b--', 'LineWidth', 1.8);
semilogx(t_data*1000, i_anion_mid_t, 'm:', 'LineWidth', 1.8);
xlabel('Time (ms)');
ylabel('Current Density (A/m^2)');
title('Current Components at Cell Midpoint');
legend('Total', 'Cation', 'Anion', 'Location', 'best');
grid on;
saveas(gcf, fullfile(output_dir, 'twoInterface_current_mid.png'));

fprintf('Plots saved to %s/\n', output_dir);

%% Display key parameters
fprintf('\n=== MODEL PARAMETERS ===\n');
fprintf('Configuration: Asymmetric (V + grounded, both with Stern layers)\n');
fprintf('Temperature: %s\n', char(model.param.get('T0')));
fprintf('Thermal voltage: %s\n', char(model.param.get('V_therm')));
fprintf('Diffusion coefficients: DA = %s, DX = %s\n', ...
    char(model.param.get('DA')), char(model.param.get('DX')));
fprintf('Bulk concentration: %s\n', char(model.param.get('cA_bulk')));
fprintf('Debye length: %s\n', char(model.param.get('xD')));
fprintf('Domain length: %s\n', char(model.param.get('L_cell')));
fprintf('Left electrode potential: %s\n', char(model.param.get('phiM_L')));
fprintf('Right electrode potential: %s (grounded)\n', char(model.param.get('phiM_R')));

%% Display key results
fprintf('\n=== KEY RESULTS (Final State) ===\n');
fprintf('Spatial profiles:\n');
fprintf('  Left electrode potential: %.4f mV\n', phi_data(1)*1000);
fprintf('  Right electrode potential: %.4f mV\n', phi_data(end)*1000);
fprintf('  Total potential drop: %.4f mV\n', (phi_data(1)-phi_data(end))*1000);
fprintf('  Center potential: %.4f mV\n', phi_data(round(end/2))*1000);
fprintf('\nConcentrations:\n');
fprintf('  Left electrode: cA = %.2f, cX = %.2f mol/m^3\n', cA_data(1), cX_data(1));
fprintf('  Right electrode: cA = %.2f, cX = %.2f mol/m^3\n', cA_data(end), cX_data(end));
fprintf('  Center: cA = %.2f, cX = %.2f mol/m^3\n', ...
    cA_data(round(end/2)), cX_data(round(end/2)));
fprintf('\nCharge conservation:\n');
fprintf('  Final Q_left: %.4e C/m^2\n', rho_L_t(end));
fprintf('  Final Q_right: %.4e C/m^2\n', rho_R_t(end));
fprintf('  Final Q_total: %.4e C/m^2 (should be ~0)\n', Q_total_t(end));
fprintf('  Charge balance error: %.2f%%\n', abs(Q_total_t(end)/(rho_L_t(end)+eps))*100);

fprintf('\n=== DONE ===\n');
fprintf('All results saved in ./%s/ directory\n', output_dir);
