% sinusoidal_one_v1.m
% Single-electrode diffuse double layer model - Sinusoidal Voltage Excitation
% Based on oneInterface_v1 with sinusoidal voltage at fixed frequency
% Parameters aligned with Mei et al. 2018 J. Phys. Chem. C 122, 194-206
% Single electrode at x=0 with Stern layer: V_dc + V_ac*sin(2*pi*f*t)
% Created: 2025-10-20

clc; clear;
close all;

%% Setup COMSOL
import com.comsol.model.*
import com.comsol.model.util.*

fprintf('Building SINGLE-ELECTRODE SINUSOIDAL model...\n');
fprintf('Configuration: One electrode at x=0 with Stern layer and sinusoidal excitation\n\n');

% Sinusoidal excitation parameters
V_dc = 0.0;               % DC bias voltage (V)
V_ac = 0.010;             % AC amplitude (10 mV)
freq = 500;               % Frequency (Hz)
omega = 2*pi*freq;        % Angular frequency (rad/s)
period = 1/freq;          % Period (s)

% Simulation time parameters
num_cycles = 10;          % Number of cycles to simulate
t_total = num_cycles * period;
pts_per_cycle = 50;       % Time points per cycle for output

fprintf('Sinusoidal Excitation Parameters:\n');
fprintf('  DC bias: %.3f V\n', V_dc);
fprintf('  AC amplitude: %.3f V (%.1f mV)\n', V_ac, V_ac*1000);
fprintf('  Frequency: %.0f Hz\n', freq);
fprintf('  Period: %.4e s\n', period);
fprintf('  Number of cycles: %d\n', num_cycles);
fprintf('  Total simulation time: %.3f s\n\n', t_total);

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

% Create study with time-dependent step only (no steady state for sinusoidal)
model.study.create('std1');
model.study('std1').create('time', 'Transient');

% Set parameters - SINGLE ELECTRODE WITH SINUSOIDAL EXCITATION
% Parameters aligned with Mei et al. 2018 JPCC paper
model.param.set('T0', '298[K]', 'Temperature (298 K from paper)');
model.param.set('V_therm', 'R_const*T0/F_const', 'Thermal voltage');
model.param.set('DA', '2e-13[m^2/s]', 'Diffusion coefficient, cation (from paper Table 2)');
model.param.set('DX', 'DA', 'Diffusion coefficient, anion (symmetric)');
model.param.set('cA_bulk', '1[mol/m^3]', 'Bulk cation concentration (0.001 M from paper)');
model.param.set('cX_bulk', 'cA_bulk', 'Bulk anion concentration');
model.param.set('zA', '+1', 'Cation charge');
model.param.set('zX', '-1', 'Anion charge');
model.param.set('Istr_bulk', '0.5*((zA^2+zX^2)*cA_bulk)', 'Bulk ionic strength');
model.param.set('eps_r', '64.4', 'Relative permittivity (propylene carbonate from paper)');
model.param.set('xD', 'sqrt(epsilon0_const*eps_r*V_therm/(2*F_const*Istr_bulk))', 'Debye length');
model.param.set('a_ion', '0.66[nm]', 'Effective ion diameter (from paper)');
model.param.set('xS', 'a_ion/2', 'Stern layer thickness (H = a/2 from paper)');
model.param.set('L_cell', '160[nm]', 'Cell length (electrolyte domain, fixed)');
model.param.set('h_max', 'L_cell/40', 'Maximum mesh element size');
model.param.set('h_max_surf', 'xD/100', 'Maximum mesh element size (electrode)');
model.param.set('Cd_GCS', 'epsilon0_const/(xD/eps_r+xS/eps_r)', 'Capacitance per unit area (GCS theory)');

% Sinusoidal voltage parameters
model.param.set('V_dc', sprintf('%.6f[V]', V_dc), 'DC bias voltage');
model.param.set('V_ac', sprintf('%.6f[V]', V_ac), 'AC voltage amplitude');
model.param.set('omega', sprintf('%.6f[rad/s]', omega), 'Angular frequency');
model.param.set('phiM', 'V_dc + V_ac*sin(omega*t)', 'Electrode potential (sinusoidal)');

% Create variables for electrode (single electrode with Stern layer)
model.component('comp1').variable.create('var1');
% Electrode variables
model.component('comp1').variable('var1').set('deltaphi', 'phiM-phi', 'Electrode-OHP potential difference');
model.component('comp1').variable('var1').set('rho_surf', 'epsilon0_const*eps_r*deltaphi/xS', 'Surface charge density');
model.component('comp1').variable('var1').set('i_charging', 'd(rho_surf,t)', 'Charging current density (A/m^2)');
% Flux and current variables
% Use tds.tflux (total flux) directly - it's the correct volume variable for bulk evaluation
model.component('comp1').variable('var1').set('i_cation', 'F_const*zA*tds.tflux_cAx', 'Cation current density (A/m^2)');
model.component('comp1').variable('var1').set('i_anion', 'F_const*zX*tds.tflux_cXx', 'Anion current density (A/m^2)');
model.component('comp1').variable('var1').set('i_total', 'i_cation+i_anion', 'Total ionic current density (A/m^2)');
% Applied voltage variable for easy extraction
model.component('comp1').variable('var1').set('V_applied', 'phiM', 'Applied voltage at electrode');

% Create geometry (1D interval from 0 to L_cell)
model.component('comp1').geom('geom1').create('i1', 'Interval');
model.component('comp1').geom('geom1').feature('i1').setIndex('coord', 'L_cell', 1);
model.component('comp1').geom('geom1').run;

% Set physics properties - Electrostatics
model.component('comp1').physics('es').feature('ccn1').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('es').feature('ccn1').set('epsilonr', {'eps_r' '0' '0' '0' 'eps_r' '0' '0' '0' 'eps_r'});

% Boundary conditions - Electrostatics
% Electrode at x=0 (boundary 1): Surface charge from Stern layer
model.component('comp1').physics('es').create('sfcd1', 'SurfaceChargeDensity', 0);
model.component('comp1').physics('es').feature('sfcd1').selection.set([1]);
model.component('comp1').physics('es').feature('sfcd1').set('rhoqs', 'rho_surf');

% Electrolyte bulk at x=L (boundary 2): Zero charge (open boundary)
model.component('comp1').physics('es').create('gnd1', 'Ground', 0);
model.component('comp1').physics('es').feature('gnd1').selection.set([2]);

% Set physics properties - Transport of Diluted Species
model.component('comp1').physics('tds').feature('sp1').setIndex('z', 'zA', 0);
model.component('comp1').physics('tds').feature('sp1').setIndex('z', 'zX', 1);
model.component('comp1').physics('tds').feature('cdm1').set('D_cA', {'DA' '0' '0' '0' 'DA' '0' '0' '0' 'DA'});
model.component('comp1').physics('tds').feature('cdm1').set('D_cX', {'DX' '0' '0' '0' 'DX' '0' '0' '0' 'DX'});
model.component('comp1').physics('tds').feature('init1').setIndex('initc', 'cA_bulk', 0);
model.component('comp1').physics('tds').feature('init1').setIndex('initc', 'cX_bulk', 1);

% Set common properties
model.common('cminpt').set('modified', {'temperature' 'T0'});

% Create mesh with refinement at electrode
model.component('comp1').mesh('mesh1').create('edg1', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg1').create('size1', 'Size');    % General mesh
model.component('comp1').mesh('mesh1').feature('edg1').create('size_surf', 'Size'); % Electrode surface

% General mesh
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmax', 'h_max');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmaxactive', true);

% Electrode refinement (boundary 1, x=0)
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_surf').selection.geom('geom1', 0);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_surf').selection.set([1]);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_surf').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_surf').set('hmax', 'h_max_surf');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_surf').set('hmaxactive', true);

model.component('comp1').mesh('mesh1').run;

% Configure time-dependent study
% Generate time list with enough points per cycle
num_time_points = num_cycles * pts_per_cycle;
t_list = linspace(0, t_total, num_time_points + 1);

fprintf('Time-dependent study configuration:\n');
fprintf('  Time range: 0 to %.4f s\n', t_total);
fprintf('  Number of output points: %d\n', length(t_list));
fprintf('  Time step: %.4e s\n', t_list(2)-t_list(1));

% Set time list
t_list_str = sprintf('%.10e ', t_list);
model.study('std1').feature('time').set('tlist', t_list_str);

fprintf('\nSingle-electrode sinusoidal model built successfully!\n');
fprintf('  Model tag: %s\n', char(model.tag()));
fprintf('  Configuration: Single electrode at x=0 with Stern layer\n');
fprintf('  Electrode: %.1f mV * sin(2*pi*%.0f*t) (with Stern layer)\n', V_ac*1000, freq);
fprintf('  Domain length: 160 nm (fixed)\n\n');

%% Run the time-dependent study
fprintf('Running time-dependent simulation...\n');
fprintf('Simulating %d cycles at %.0f Hz...\n', num_cycles, freq);
tic;
model.study('std1').run();
elapsed_time = toc;
fprintf('Simulation completed in %.2f seconds\n', elapsed_time);

%% Extract results along the domain (multiple time snapshots)
fprintf('\nExtracting spatial profiles at selected time points...\n');

% Select a few time snapshots to save (first, middle, last, and a few in between)
snapshot_indices = unique(round(linspace(1, length(t_list), min(10, length(t_list)))));
num_snapshots = length(snapshot_indices);

fprintf('  Saving %d spatial snapshots\n', num_snapshots);

% Extract first snapshot to get spatial dimensions
eval_result = mpheval(model, {'x', 'phi', 'cA', 'cX', 'es.normE', 'es.rhoq'}, ...
    'edim', 1, 'solnum', snapshot_indices(1));
num_x_points = length(eval_result.d1);

% Preallocate arrays
x_data = eval_result.d1;
phi_snapshots = zeros(num_x_points, num_snapshots);
cA_snapshots = zeros(num_x_points, num_snapshots);
cX_snapshots = zeros(num_x_points, num_snapshots);
Es_snapshots = zeros(num_x_points, num_snapshots);
rho_snapshots = zeros(num_x_points, num_snapshots);
t_snapshots = t_list(snapshot_indices);

% Extract all snapshots
for i = 1:num_snapshots
    result = mpheval(model, {'phi', 'cA', 'cX', 'es.normE', 'es.rhoq'}, ...
        'edim', 1, 'solnum', snapshot_indices(i));
    phi_snapshots(:, i) = result.d1;
    cA_snapshots(:, i) = result.d2;
    cX_snapshots(:, i) = result.d3;
    Es_snapshots(:, i) = result.d4;
    rho_snapshots(:, i) = result.d5;
end

fprintf('  Number of spatial points: %d\n', num_x_points);

%% Extract time series data at BOTH electrodes
fprintf('\nExtracting time evolution at left electrode (x=0)...\n');
time_L = mpheval(model, {'t', 'V_applied', 'phi', 'cA', 'cX', 'i_charging_L', 'i_cation', 'i_anion', 'i_total', 'rho_surf_L'}, ...
    'edim', 0, 'selection', 1);
t_data = time_L.d1;
V_applied_t = time_L.d2;
phi_L_t = time_L.d3;
cA_L_t = time_L.d4;
cX_L_t = time_L.d5;
i_charging_L_t = time_L.d6;
i_cation_L_t = time_L.d7;
i_anion_L_t = time_L.d8;
i_total_L_t = time_L.d9;
rho_L_t = time_L.d10;

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

% Get x coordinates from first solution to find midpoint index
eval_x = mpheval(model, 'x', 'edim', 1, 'solnum', 1);
x_coords = eval_x.d1;
[~, mid_idx] = min(abs(x_coords - x_mid)); % Find closest point to midpoint

fprintf('  Using mesh point at x = %.4e m (index %d)\n', x_coords(mid_idx), mid_idx);

% Preallocate arrays
num_times = length(t_data);
i_total_mid_t = zeros(num_times, 1);
i_cation_mid_t = zeros(num_times, 1);
i_anion_mid_t = zeros(num_times, 1);

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
output_dir = 'rst/sinusoidal';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('\nCreated output directory: %s\n', output_dir);
end

%% Save COMSOL model file
model_file = fullfile(output_dir, 'sinusoidal.mph');
mphsave(model, model_file);
fprintf('COMSOL model saved to: %s\n', model_file);

%% Save results to file
results_file = fullfile(output_dir, 'sinusoidal_results.mat');
save(results_file, 'x_data', 'phi_snapshots', 'cA_snapshots', 'cX_snapshots', ...
    'Es_snapshots', 'rho_snapshots', 't_snapshots', ...
    't_data', 'V_applied_t', 'phi_L_t', 'cA_L_t', 'cX_L_t', 'phi_R_t', 'cA_R_t', 'cX_R_t', ...
    'i_charging_L_t', 'i_cation_L_t', 'i_anion_L_t', 'i_total_L_t', 'rho_L_t', ...
    'i_charging_R_t', 'i_cation_R_t', 'i_anion_R_t', 'i_total_R_t', 'rho_R_t', ...
    'x_mid', 'i_cation_mid_t', 'i_anion_mid_t', 'i_total_mid_t', ...
    'Q_total_t', 'V_dc', 'V_ac', 'freq', 'omega', 'period', 'num_cycles', 'elapsed_time');
fprintf('Results saved to: %s\n', results_file);

%% Save time evolution data
csv_time_L = fullfile(output_dir, 'sinusoidal_time_left.csv');
time_L_table = table(t_data*1000, V_applied_t*1000, phi_L_t*1000, cA_L_t, cX_L_t, ...
    i_charging_L_t, i_cation_L_t, i_anion_L_t, i_total_L_t, rho_L_t, ...
    'VariableNames', {'t_ms', 'V_applied_mV', 'phi_mV', 'cA_molm3', 'cX_molm3', ...
    'i_charging_Am2', 'i_cation_Am2', 'i_anion_Am2', 'i_total_Am2', 'rho_surf_Cm2'});
writetable(time_L_table, csv_time_L);
fprintf('Left electrode time data saved to: %s\n', csv_time_L);

csv_time_R = fullfile(output_dir, 'sinusoidal_time_right.csv');
time_R_table = table(t_data*1000, phi_R_t*1000, cA_R_t, cX_R_t, ...
    i_charging_R_t, i_cation_R_t, i_anion_R_t, i_total_R_t, rho_R_t, ...
    'VariableNames', {'t_ms', 'phi_mV', 'cA_molm3', 'cX_molm3', ...
    'i_charging_Am2', 'i_cation_Am2', 'i_anion_Am2', 'i_total_Am2', 'rho_surf_Cm2'});
writetable(time_R_table, csv_time_R);
fprintf('Right electrode time data saved to: %s\n', csv_time_R);

csv_time_mid = fullfile(output_dir, 'sinusoidal_time_mid.csv');
time_mid_table = table(t_data*1000, V_applied_t*1000, i_total_mid_t, i_cation_mid_t, i_anion_mid_t, ...
    'VariableNames', {'t_ms', 'V_applied_mV', 'i_total_Am2', 'i_cation_Am2', 'i_anion_Am2'});
writetable(time_mid_table, csv_time_mid);
fprintf('Midpoint time data saved to: %s\n', csv_time_mid);

%% ========================================================================
%% IMPEDANCE ANALYSIS
%% ========================================================================

fprintf('\n=== IMPEDANCE ANALYSIS ===\n');
fprintf('Excitation frequency: %.0f Hz\n', freq);
fprintf('Excitation amplitude: %.3f V (%.1f mV)\n', V_ac, V_ac*1000);
fprintf('Angular frequency: %.3f rad/s\n', omega);
fprintf('Period: %.4e s\n\n', period);

%% Extract last 2 cycles for steady periodic analysis
% Ignore transient startup behavior
idx_last_cycles = t_data >= (max(t_data) - 2*period);
t_ss = t_data(idx_last_cycles);
V_ss = V_applied_t(idx_last_cycles);
I_ss = i_total_mid_t(idx_last_cycles);

num_points_ss = length(t_ss);
fprintf('Analyzing last 2 cycles:\n');
fprintf('  Number of points: %d\n', num_points_ss);
fprintf('  Time range: %.4e to %.4e s\n', min(t_ss), max(t_ss));

%% Calculate impedance using FFT method
% Voltage amplitude (should be close to V_ac)
V_amplitude = (max(V_ss) - min(V_ss)) / 2;
I_amplitude = (max(I_ss) - min(I_ss)) / 2;

fprintf('\nVoltage amplitude: %.4e V (%.3f mV)\n', V_amplitude, V_amplitude*1000);
fprintf('Current amplitude: %.4e A/m^2\n', I_amplitude);

% Center the signals (remove DC offset)
V_mean = mean(V_ss);
I_mean = mean(I_ss);
V_centered = V_ss - V_mean;
I_centered = I_ss - I_mean;

% Apply Hanning window to reduce spectral leakage
N = length(V_centered);
n = (0:N-1)';
window = 0.5 * (1 - cos(2*pi*n/(N-1)));
V_windowed = V_centered .* window;
I_windowed = I_centered .* window;

% FFT
NFFT = 2^nextpow2(length(V_windowed));
V_fft = fft(V_windowed, NFFT);
I_fft = fft(I_windowed, NFFT);

% Frequency vector
Fs = 1 / (t_ss(2) - t_ss(1));  % Sampling frequency
f_vector = Fs * (0:(NFFT/2-1)) / NFFT;

% Find the peak at excitation frequency
[~, idx_peak] = min(abs(f_vector - freq));

% Complex impedance at excitation frequency
Z_complex = V_fft(idx_peak) / I_fft(idx_peak);
Z_mag = abs(Z_complex);
phase_rad = -angle(Z_complex);  % Negative because Z = V/I
phase_deg = phase_rad * 180/pi;

% Real and imaginary parts
Z_real = real(Z_complex);
Z_imag = imag(Z_complex);

fprintf('\nImpedance Results (FFT Analysis at %.0f Hz):\n', f_vector(idx_peak));
fprintf('  Complex impedance: Z = %.4e %+.4ej Ohm*m^2\n', Z_real, Z_imag);
fprintf('  Magnitude: |Z| = %.4e Ohm*m^2\n', Z_mag);
fprintf('  Phase angle: %.3f rad (%.2f degrees)\n', phase_rad, phase_deg);

%% Physical interpretation
fprintf('\n=== PHYSICAL INTERPRETATION ===\n');

if abs(phase_deg) < 10
    fprintf('Phase ~ 0 deg: Resistive behavior dominates\n');
elseif phase_deg > 80
    fprintf('Phase ~ 90 deg: Capacitive behavior dominates\n');
    C_eff = -1 / (omega * Z_imag);
    fprintf('  Effective capacitance: %.4e F/m^2\n', C_eff);
elseif phase_deg < -80
    fprintf('Phase ~ -90 deg: Inductive behavior dominates\n');
else
    fprintf('Phase ~ %.1f deg: Mixed resistive-capacitive behavior\n', phase_deg);
    if Z_imag < 0
        C_eff = -1 / (omega * Z_imag);
        fprintf('  Resistive component: %.4e Ohm*m^2\n', Z_real);
        fprintf('  Capacitive component: %.4e F/m^2\n', C_eff);
    end
end

%% Compare with theoretical capacitance
if Z_imag < 0
    C_measured = -1 / (omega * Z_imag);
    fprintf('\n=== CAPACITANCE COMPARISON ===\n');
    fprintf('Measured capacitance: %.4e F/m^2\n', C_measured);

    % Calculate theoretical Stern layer capacitance
    epsilon0 = 8.854e-12;  % F/m
    eps_S = 10;            % Relative permittivity of Stern layer
    xS = 0.2e-9;           % Stern layer thickness (m)
    C_stern = epsilon0 * eps_S / xS;

    fprintf('Theoretical Stern capacitance: %.4e F/m^2\n', C_stern);
    fprintf('Ratio (measured/theory): %.3f\n', C_measured / C_stern);
end

%% Save impedance results
impedance_results = struct();
impedance_results.freq = freq;
impedance_results.omega = omega;
impedance_results.V_ac = V_ac;
impedance_results.V_amplitude = V_amplitude;
impedance_results.I_amplitude = I_amplitude;
impedance_results.Z_magnitude = Z_mag;
impedance_results.phase_rad = phase_rad;
impedance_results.phase_deg = phase_deg;
impedance_results.Z_real = Z_real;
impedance_results.Z_imag = Z_imag;
impedance_results.Z_complex = Z_complex;

save(fullfile(output_dir, 'impedance_results.mat'), 'impedance_results');
fprintf('Impedance results saved to %s/impedance_results.mat\n', output_dir);

%% ========================================================================
%% PLOTTING
%% ========================================================================

fprintf('\n=== GENERATING PLOTS ===\n');

% Figure 1: Time evolution of applied voltage and response (full simulation)
figure(1);
subplot(3,1,1);
plot(t_data*1000, V_applied_t*1000, 'k-', 'LineWidth', 2); hold on;
plot(t_data*1000, phi_L_t*1000, 'r-', 'LineWidth', 1.5);
plot(t_data*1000, phi_R_t*1000, 'b-', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Potential (mV)');
title(sprintf('Voltage and Electrode Potentials (%.0f Hz)', freq));
legend('Applied V', 'Left \phi', 'Right \phi', 'Location', 'best');
grid on;

subplot(3,1,2);
plot(t_data*1000, i_total_mid_t, 'k-', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Current Density (A/m^2)');
title('Total Current at Cell Midpoint');
grid on;

subplot(3,1,3);
plot(t_data*1000, rho_L_t, 'r-', 'LineWidth', 2); hold on;
plot(t_data*1000, rho_R_t, 'b-', 'LineWidth', 2);
plot(t_data*1000, Q_total_t, 'k--', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Surface Charge Density (C/m^2)');
title('Surface Charge Density');
legend('Left', 'Right', 'Total', 'Location', 'best');
grid on;

saveas(gcf, fullfile(output_dir, 'sinusoidal_time_evolution.png'));

% Figure 2: Last 2 cycles - Voltage and Current waveforms
figure(2);
subplot(2,1,1);
plot(t_ss*1000, V_ss*1000, 'r-', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title(sprintf('Applied Voltage - Last 2 Cycles (%.0f Hz)', freq));
grid on;

subplot(2,1,2);
plot(t_ss*1000, I_ss, 'b-', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Current (A/m^2)');
title('Midpoint Current Response - Last 2 Cycles');
grid on;

saveas(gcf, fullfile(output_dir, 'impedance_waveforms.png'));

% Figure 3: V-I Lissajous curve
figure(3);
plot(V_ss*1000, I_ss, 'b-', 'LineWidth', 2);
xlabel('Voltage (mV)');
ylabel('Current (A/m^2)');
title(sprintf('V-I Lissajous (Z = %.2e Ohm*m^2, \\phi = %.1f deg)', Z_mag, phase_deg));
grid on;
xlim([-10 10]);
ylim([-0.15 0.15]);
saveas(gcf, fullfile(output_dir, 'impedance_lissajous.png'));

% Figure 4: Nyquist plot (single frequency point)
figure(4);
plot(Z_real, -Z_imag, 'ro', 'MarkerSize', 12, 'LineWidth', 2);
xlabel('Z'' (Ohm*m^2)');
ylabel('-Z" (Ohm*m^2)');
title(sprintf('Nyquist Plot at %.0f Hz', freq));
grid on;

% Set axis limits with some padding around the data point
Z_max = max(abs([Z_real, Z_imag]));
if Z_max > 0
    margin = Z_max * 0.2;  % 20% margin
    xlim([min(0, Z_real) - margin, Z_real + margin]);
    ylim([min(0, -Z_imag) - margin, -Z_imag + margin]);
end
axis equal;

saveas(gcf, fullfile(output_dir, 'impedance_nyquist.png'));

% Figure 5: Bode plot (single frequency point)
figure(5);
subplot(2,1,1);
semilogx(freq, Z_mag, 'bo', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('|Z| (Ohm*m^2)');
title('Bode Plot - Magnitude');
grid on;
xlim([freq/10 freq*10]);
ylim([0 Z_mag*1.5]);

subplot(2,1,2);
semilogx(freq, phase_deg, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
title('Bode Plot - Phase');
grid on;
xlim([freq/10 freq*10]);
ylim([-5 100]);

saveas(gcf, fullfile(output_dir, 'impedance_bode.png'));

% Figure 6: Spatial profiles at different times
figure(6);
subplot(2,2,1);
hold on;
for i = 1:num_snapshots
    plot(x_data*1e9, phi_snapshots(:,i)*1000, 'LineWidth', 1.5);
end
xlabel('Position (nm)');
ylabel('Potential (mV)');
title('Electric Potential Snapshots');
colorbar;
grid on;

subplot(2,2,2);
hold on;
for i = 1:num_snapshots
    plot(x_data*1e9, cA_snapshots(:,i), 'r-', 'LineWidth', 1);
    plot(x_data*1e9, cX_snapshots(:,i), 'b-', 'LineWidth', 1);
end
xlabel('Position (nm)');
ylabel('Concentration (mol/m^3)');
title('Ion Concentration Snapshots');
grid on;

subplot(2,2,3);
hold on;
for i = 1:num_snapshots
    plot(x_data*1e9, Es_snapshots(:,i)/1e6, 'LineWidth', 1);
end
xlabel('Position (nm)');
ylabel('E-field (MV/m)');
title('Electric Field Snapshots');
grid on;

subplot(2,2,4);
hold on;
for i = 1:num_snapshots
    plot(x_data*1e9, rho_snapshots(:,i), 'LineWidth', 1);
end
xlabel('Position (nm)');
ylabel('Charge density (C/m^3)');
title('Space Charge Density Snapshots');
grid on;

saveas(gcf, fullfile(output_dir, 'sinusoidal_spatial_snapshots.png'));

fprintf('All plots saved to %s/\n', output_dir);

%% Display key parameters
fprintf('\n=== MODEL PARAMETERS ===\n');
fprintf('Configuration: Asymmetric sinusoidal (V_ac*sin(wt) + grounded)\n');
fprintf('Temperature: %s\n', char(model.param.get('T0')));
fprintf('Thermal voltage: %s\n', char(model.param.get('V_therm')));
fprintf('Diffusion coefficients: DA = %s, DX = %s\n', ...
    char(model.param.get('DA')), char(model.param.get('DX')));
fprintf('Bulk concentration: %s\n', char(model.param.get('cA_bulk')));
fprintf('Debye length: %s\n', char(model.param.get('xD')));
fprintf('Domain length: %s\n', char(model.param.get('L_cell')));
fprintf('Sinusoidal parameters:\n');
fprintf('  DC bias: %.3f V\n', V_dc);
fprintf('  AC amplitude: %.3f V (%.1f mV)\n', V_ac, V_ac*1000);
fprintf('  Frequency: %.0f Hz\n', freq);
fprintf('  Number of cycles: %d\n', num_cycles);

%% Display key results
fprintf('\n=== KEY RESULTS ===\n');
fprintf('Applied voltage range: %.4f to %.4f mV\n', min(V_applied_t)*1000, max(V_applied_t)*1000);
fprintf('Midpoint current range: %.4e to %.4e A/m^2\n', min(i_total_mid_t), max(i_total_mid_t));
fprintf('Charging current range (left): %.4e to %.4e A/m^2\n', min(i_charging_L_t), max(i_charging_L_t));
fprintf('Surface charge range (left): %.4e to %.4e C/m^2\n', min(rho_L_t), max(rho_L_t));
fprintf('Surface charge range (right): %.4e to %.4e C/m^2\n', min(rho_R_t), max(rho_R_t));
fprintf('\nCharge conservation:\n');
fprintf('  Q_total range: %.4e to %.4e C/m^2\n', min(Q_total_t), max(Q_total_t));
fprintf('  Max charge imbalance: %.4e C/m^2 (should be ~0)\n', max(abs(Q_total_t)));

fprintf('\n=== DONE ===\n');
fprintf('All results saved in ./%s/ directory\n', output_dir);
