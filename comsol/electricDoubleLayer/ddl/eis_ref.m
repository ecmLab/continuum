% eis_ref.m
% Electrochemical Impedance Spectroscopy (EIS) - Reference Paper Implementation
% Based on: Mei et al., J. Phys. Chem. C 2018, 122, 194-206
% "Physical Interpretations of Nyquist Plots for EDLC Electrodes and Devices"
%
% Implementation follows Supporting Information equations S.1-S.15
% Uses least-squares fitting approach (not FFT) as in paper
% Single electrode setup (blocking electrode at x=0)
%
% PARAMETERS (following paper Table 1 and baseline case 5 from Table 2):
%   - Temperature: T = 298 K (25°C)
%   - Electrolyte: ε_r = 64.4 (propylene carbonate)
%   - Ion valency: z = 1 (binary symmetric ±1)
%   - Ion diameter: a = 0.66 nm → Stern layer: H = a/2 = 0.33 nm
%   - Diffusion coefficient: D = 2×10⁻¹³ m²/s (both ions, symmetric)
%   - Bulk concentration: c∞ = 1 mol/m³ (0.001 mol/L)
%   - Domain length: L = 160 nm
%   - DC bias: ψ_dc = 0.0 V (EIS at equilibrium)
%   - AC amplitude: ψ₀ = 5 mV
%   - Frequency range: 0.1 Hz to 5 MHz
%
% NOTE: Electrode conductivity and thickness are NOT included (electrolyte-only model)
% NOTE: V_dc = 0 eliminates DC transient (uniform initial condition = equilibrium)
%
% Created: 2025-10-20
% Updated: 2025-10-20 - Parameters matched to paper Table 1
% Updated: 2025-10-20 - Fixed V_dc = 0 and syntax errors

clc; clear;
close all;

%% EIS Configuration (following paper methodology)
fprintf('=== EIS REFERENCE IMPLEMENTATION ===\n');
fprintf('Based on: Mei et al., J. Phys. Chem. C 2018, 122, 194-206\n\n');

% Frequency range: 0.1 Hz to 5 MHz (matching paper Table 1)
% Paper states: f = 0.1 - 5×10^6 Hz
frequencies = logspace(-1, log10(5e6), 25);  % 0.1 Hz to 5 MHz, 25 points
num_freqs = length(frequencies);

% Excitation parameters
% Note: V_dc = 0 for EIS (standard practice - measure around equilibrium)
% Paper's V_dc = 0.3 V is likely for galvanostatic cycling, not EIS
V_dc = 0.0;               % DC bias voltage (V) - EIS at equilibrium
V_ac = 0.005;             % AC amplitude (5 mV) - paper uses 5 mV
num_cycles = 10;          % Number of cycles per frequency
pts_per_cycle = 50;       % Time points per cycle (for least-squares fitting)

fprintf('EIS Configuration:\n');
fprintf('  Frequency range: %.1f Hz to %.0e Hz (%d points, log-spaced)\n', ...
    frequencies(1), frequencies(end), num_freqs);
fprintf('  DC bias: %.3f V\n', V_dc);
fprintf('  AC amplitude: %.3f V (%.1f mV)\n', V_ac, V_ac*1000);
fprintf('  Cycles per frequency: %d\n', num_cycles);
fprintf('  Time points per cycle: %d\n\n', pts_per_cycle);

%% Create output directory
output_dir = 'rst/eisRef';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n\n', output_dir);
end

%% Initialize storage for EIS results
Z_real_all = zeros(num_freqs, 1);
Z_imag_all = zeros(num_freqs, 1);
Z_mag_all = zeros(num_freqs, 1);
phase_deg_all = zeros(num_freqs, 1);

% Also store fitted current components for validation
I_dc_all = zeros(num_freqs, 1);
I_sin_all = zeros(num_freqs, 1);
I_cos_all = zeros(num_freqs, 1);
fit_R2_all = zeros(num_freqs, 1);

%% Setup COMSOL - Single Electrode Configuration
import com.comsol.model.*
import com.comsol.model.util.*

fprintf('Building COMSOL model (single electrode)...\n');

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

% Create study with time-dependent step
model.study.create('std1');
model.study('std1').create('time', 'Transient');

% Set parameters - SINGLE ELECTRODE CONFIGURATION (paper eq. S.5a, S.7-S.10)
% Following paper Table 1 and baseline case from Table 2 (case 5)
model.param.set('T0', '25[degC]', 'Temperature (T=298K)');
model.param.set('V_therm', 'R_const*T0/F_const', 'Thermal voltage');

% Electrolyte properties (paper Table 1)
model.param.set('eps_r', '64.4', 'Relative permittivity (propylene carbonate, PC)');
model.param.set('z', '1', 'Ion valency (binary symmetric ±1)');

% Ion diameter and Stern layer (paper approach: H = a/2)
model.param.set('a', '0.66[nm]', 'Ion diameter (baseline from Table 2)');
model.param.set('H', 'a/2', 'Stern layer thickness (H = a/2, paper assumption)');

% Diffusion coefficients (paper Table 1 baseline, case 5)
% Paper uses D = 2e-13 m²/s as baseline for symmetric ions
model.param.set('D', '2e-13[m^2/s]', 'Diffusion coefficient (both ions, symmetric)');
model.param.set('DA', 'D', 'Diffusion coefficient, cation');
model.param.set('DX', 'D', 'Diffusion coefficient, anion');

% Concentrations - binary symmetric electrolyte (paper baseline case 5)
model.param.set('c_inf', '1[mol/m^3]', 'Bulk concentration (0.001 mol/L from Table 2)');
model.param.set('cA_bulk', 'c_inf', 'Bulk cation concentration');
model.param.set('cX_bulk', 'c_inf', 'Bulk anion concentration');

% Ion charges (symmetric ±z)
model.param.set('zA', '+z', 'Cation charge');
model.param.set('zX', '-z', 'Anion charge');

% Derived parameters
model.param.set('Istr_bulk', '0.5*((zA^2+zX^2)*c_inf)', 'Bulk ionic strength');
model.param.set('lambda_D', 'sqrt(epsilon0_const*eps_r*V_therm/(2*F_const*Istr_bulk))', 'Debye length');

% Domain length (paper Table 1: L ranges 40-1600 nm, baseline case 5: 160 nm)
model.param.set('L', '160[nm]', 'Electrolyte domain thickness (paper baseline)');

% Mesh parameters
model.param.set('h_max', 'L/40', 'Maximum mesh element size');
model.param.set('h_max_surf', 'lambda_D/100', 'Maximum mesh element size (electrode)');

% Stern layer capacitance (paper eq. S.16-S.17)
% Note: Paper doesn't explicitly state eps_S, commonly assumed same as electrolyte or ~10
model.param.set('eps_S', '10', 'Relative permittivity of Stern layer');
model.param.set('Cd_GCS', 'epsilon0_const/(lambda_D/eps_r+H/eps_S)', 'Capacitance per unit area (GCS theory)');

% Sinusoidal voltage at electrode (paper eq. S.5a)
model.param.set('V_dc', sprintf('%.6f[V]', V_dc), 'DC bias voltage');
model.param.set('V_ac', sprintf('%.6f[V]', V_ac), 'AC voltage amplitude');
model.param.set('omega', '2*pi*1e2[rad/s]', 'Angular frequency (placeholder)');
model.param.set('phiM', 'V_dc + V_ac*sin(omega*t)', 'Electrode potential (sinusoidal)');

% Create variables for electrode-electrolyte interface
model.component('comp1').variable.create('var1');

% Electrode potential difference and surface charge (paper eq. S.1, S.4)
model.component('comp1').variable('var1').set('deltaphi', 'phiM-phi', 'Electrode-OHP potential difference');
model.component('comp1').variable('var1').set('rho_surf', 'epsilon0_const*eps_S*deltaphi/H', 'Surface charge density (Stern layer)');
model.component('comp1').variable('var1').set('i_charging', 'd(rho_surf,t)', 'Charging current density (A/m^2)');

% Ionic current densities (paper eq. S.3)
model.component('comp1').variable('var1').set('i_cation', 'F_const*zA*tds.tflux_cAx', 'Cation current density (A/m^2)');
model.component('comp1').variable('var1').set('i_anion', 'F_const*zX*tds.tflux_cXx', 'Anion current density (A/m^2)');
model.component('comp1').variable('var1').set('i_total', 'i_cation+i_anion', 'Total ionic current density (A/m^2)');

% Applied voltage (for impedance calculation)
model.component('comp1').variable('var1').set('V_applied', 'phiM', 'Applied voltage');

% Create geometry (1D interval from 0 to L)
model.component('comp1').geom('geom1').create('i1', 'Interval');
model.component('comp1').geom('geom1').feature('i1').setIndex('coord', 'L', 1);
model.component('comp1').geom('geom1').run;

% Set physics properties - Electrostatics
model.component('comp1').physics('es').feature('ccn1').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('es').feature('ccn1').set('epsilonr', {'eps_r' '0' '0' '0' 'eps_r' '0' '0' '0' 'eps_r'});

% Boundary conditions - Electrostatics
% LEFT (x=0): Blocking electrode with surface charge (paper eq. S.7)
model.component('comp1').physics('es').create('sfcd1', 'SurfaceChargeDensity', 0);
model.component('comp1').physics('es').feature('sfcd1').selection.set([1]);
model.component('comp1').physics('es').feature('sfcd1').set('rhoqs', 'rho_surf');

% RIGHT (x=L): Bulk electrolyte, zero potential reference (paper eq. S.8)
model.component('comp1').physics('es').create('gnd1', 'Ground', 0);
model.component('comp1').physics('es').feature('gnd1').selection.set([2]);

% Boundary conditions - Transport of Diluted Species
% LEFT (x=0): No flux (blocking electrode, paper eq. S.9)
% (Default is no flux, no need to add BC)

% RIGHT (x=L): Bulk concentration (paper eq. S.10)
model.component('comp1').physics('tds').create('conc1', 'Concentration', 0);
model.component('comp1').physics('tds').feature('conc1').selection.set([2]);
model.component('comp1').physics('tds').feature('conc1').setIndex('species', true, 0);
model.component('comp1').physics('tds').feature('conc1').setIndex('species', true, 1);
model.component('comp1').physics('tds').feature('conc1').setIndex('c0', 'cA_bulk', 0);
model.component('comp1').physics('tds').feature('conc1').setIndex('c0', 'cX_bulk', 1);

% Set physics properties - Transport of Diluted Species
model.component('comp1').physics('tds').feature('sp1').setIndex('z', 'zA', 0);
model.component('comp1').physics('tds').feature('sp1').setIndex('z', 'zX', 1);
model.component('comp1').physics('tds').feature('cdm1').set('D_cA', {'DA' '0' '0' '0' 'DA' '0' '0' '0' 'DA'});
model.component('comp1').physics('tds').feature('cdm1').set('D_cX', {'DX' '0' '0' '0' 'DX' '0' '0' '0' 'DX'});

% Initial conditions - uniform bulk concentration
model.component('comp1').physics('tds').feature('init1').setIndex('initc', 'cA_bulk', 0);
model.component('comp1').physics('tds').feature('init1').setIndex('initc', 'cX_bulk', 1);

% Set common properties (temperature)
model.common('cminpt').set('modified', {'temperature' 'T0'});

% Create mesh with refinement at electrode (x=0)
model.component('comp1').mesh('mesh1').create('edg1', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg1').create('size_elec', 'Size');

% General mesh
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmax', 'h_max');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmaxactive', true);

% Refined mesh at electrode (point 1, x=0)
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_elec').selection.geom('geom1', 0);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_elec').selection.set([1]);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_elec').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_elec').set('hmax', 'h_max_surf');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_elec').set('hmaxactive', true);

model.component('comp1').mesh('mesh1').run;

fprintf('Base model built successfully!\n');
fprintf('  Configuration: Single blocking electrode at x=0\n');
fprintf('  Boundary conditions: Surface charge (x=0), Ground + Bulk conc (x=L)\n\n');

%% Frequency sweep loop (high to low, standard EIS convention)
fprintf('Starting frequency sweep (%d frequencies)...\n', num_freqs);
fprintf('Sweeping from HIGH to LOW frequency (standard EIS)\n\n');

for freq_idx = 1:num_freqs
    freq = frequencies(freq_idx);
    omega = 2*pi*freq;
    period = 1/freq;
    t_total = num_cycles * period;

    fprintf('========================================\n');
    fprintf('Frequency %d/%d: %.2e Hz\n', freq_idx, num_freqs, freq);
    fprintf('========================================\n');
    fprintf('  Period: %.4e s\n', period);
    fprintf('  Total time: %.4e s (%d cycles)\n', t_total, num_cycles);

    %% Update frequency parameter
    model.param.set('omega', sprintf('%.10e[rad/s]', omega), 'Angular frequency');

    %% Configure time-dependent study
    num_time_points = num_cycles * pts_per_cycle;
    t_list = linspace(0, t_total, num_time_points + 1);

    fprintf('  Time points: %d\n', length(t_list));
    fprintf('  Time step: %.4e s\n', t_list(2)-t_list(1));

    % Set time list
    t_list_str = sprintf('%.10e ', t_list);
    model.study('std1').feature('time').set('tlist', t_list_str);

    %% Run simulation
    fprintf('  Running simulation...\n');
    tic;
    model.study('std1').run();
    elapsed = toc;
    fprintf('  Simulation completed in %.2f seconds\n', elapsed);

    %% Extract time series data at ELECTRODE (x=0, boundary point 1)
    % This follows paper's approach of measuring electrode current
    fprintf('  Extracting data at electrode (x=0)...\n');

    % Get time, voltage, and current at electrode
    time_data = mpheval(model, {'t', 'V_applied', 'i_total'}, 'edim', 0, 'selection', 1);
    t_data = time_data.d1;
    V_applied_t = time_data.d2;
    I_electrode_t = time_data.d3;  % Current density at electrode (A/m^2)

    %% Impedance analysis using LEAST-SQUARES FITTING (paper method)
    % Paper: fit I(t) = I_dc + I_sin*sin(ωt) + I_cos*cos(ωt)
    % Then: Z = V_amp / (I_sin - j*I_cos)

    fprintf('  Analyzing impedance (least-squares fitting)...\n');

    % Use last 2 cycles for steady-state analysis
    idx_last_cycles = t_data >= (max(t_data) - 2*period);
    t_ss = t_data(idx_last_cycles);
    V_ss = V_applied_t(idx_last_cycles);
    I_ss = I_electrode_t(idx_last_cycles);

    % Build design matrix for least-squares
    % Model: I(t) = β0 + β1*sin(ωt) + β2*cos(ωt)
    A_matrix = [ones(length(t_ss), 1), sin(omega*t_ss), cos(omega*t_ss)];

    % Solve least-squares: β = (A'A)^(-1) A'I
    beta = A_matrix \ I_ss;
    I_dc_fit = beta(1);
    I_sin_fit = beta(2);
    I_cos_fit = beta(3);

    % Calculate goodness of fit (R^2)
    I_fit = A_matrix * beta;
    SS_res = sum((I_ss - I_fit).^2);
    SS_tot = sum((I_ss - mean(I_ss)).^2);
    R2 = 1 - SS_res/SS_tot;

    % Complex impedance (paper equation)
    % Z = V_amp / (I_sin - j*I_cos)
    Z_complex = V_ac / (I_sin_fit - 1j*I_cos_fit);

    Z_real = real(Z_complex);
    Z_imag = imag(Z_complex);
    Z_mag = abs(Z_complex);
    phase_rad = angle(Z_complex);
    phase_deg = phase_rad * 180/pi;

    % Store results
    Z_real_all(freq_idx) = Z_real;
    Z_imag_all(freq_idx) = Z_imag;
    Z_mag_all(freq_idx) = Z_mag;
    phase_deg_all(freq_idx) = phase_deg;
    I_dc_all(freq_idx) = I_dc_fit;
    I_sin_all(freq_idx) = I_sin_fit;
    I_cos_all(freq_idx) = I_cos_fit;
    fit_R2_all(freq_idx) = R2;

    fprintf('\n  Least-squares fit:\n');
    fprintf('    I_dc = %.4e A/m^2\n', I_dc_fit);
    fprintf('    I_sin = %.4e A/m^2\n', I_sin_fit);
    fprintf('    I_cos = %.4e A/m^2\n', I_cos_fit);
    fprintf('    R^2 = %.6f\n', R2);

    fprintf('\n  Impedance:\n');
    fprintf('    Z = %.4e %+.4ej Ohm*m^2\n', Z_real, Z_imag);
    fprintf('    |Z| = %.4e Ohm*m^2\n', Z_mag);
    fprintf('    Phase = %.2f deg\n', phase_deg);

    fprintf('\nProgress: %d/%d frequencies completed\n\n', freq_idx, num_freqs);
end

fprintf('========================================\n');
fprintf('=== FREQUENCY SWEEP COMPLETED ===\n');
fprintf('========================================\n\n');

%% Save EIS data
fprintf('Saving results...\n');

% Save CSV with all data including fit quality
csv_file = fullfile(output_dir, 'eis_ref_data.csv');
eis_table = table(frequencies', Z_real_all, Z_imag_all, Z_mag_all, phase_deg_all, ...
    I_dc_all, I_sin_all, I_cos_all, fit_R2_all, ...
    'VariableNames', {'freq_Hz', 'Z_real_Ohm_m2', 'Z_imag_Ohm_m2', 'Z_mag_Ohm_m2', ...
    'phase_deg', 'I_dc_A_m2', 'I_sin_A_m2', 'I_cos_A_m2', 'fit_R2'});
writetable(eis_table, csv_file);
fprintf('  CSV: %s\n', csv_file);

% Save MAT file
mat_file = fullfile(output_dir, 'eis_ref_results.mat');
save(mat_file, 'frequencies', 'Z_real_all', 'Z_imag_all', 'Z_mag_all', 'phase_deg_all', ...
    'I_dc_all', 'I_sin_all', 'I_cos_all', 'fit_R2_all', ...
    'V_dc', 'V_ac', 'num_cycles', 'pts_per_cycle');
fprintf('  MAT: %s\n', mat_file);

% Save COMSOL model
mph_file = fullfile(output_dir, 'eis_ref_solved.mph');
mphsave(model, mph_file);
fprintf('  COMSOL: %s\n', mph_file);

%% Generate plots
fprintf('\nGenerating plots...\n');

% Figure 1: Nyquist plot with feature identification
figure(1);
clf;
plot(Z_real_all, -Z_imag_all, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;

% Mark high and low frequency points
plot(Z_real_all(1), -Z_imag_all(1), 'rs', 'MarkerSize', 12, 'LineWidth', 2);
plot(Z_real_all(end), -Z_imag_all(end), 'gs', 'MarkerSize', 12, 'LineWidth', 2);

% Add frequency labels at select points
label_indices = [1, round(num_freqs/4), round(num_freqs/2), round(3*num_freqs/4), num_freqs];
for i = label_indices
    text(Z_real_all(i), -Z_imag_all(i), sprintf('  %.0e Hz', frequencies(i)), ...
        'FontSize', 8, 'VerticalAlignment', 'bottom');
end

xlabel('Z'' (Ohm*m^2)');
ylabel('-Z" (Ohm*m^2)');
title('Nyquist Plot - Reference Implementation (Mei et al. 2018)');
legend('EIS Data', sprintf('High Freq (%.0e Hz)', frequencies(1)), ...
    sprintf('Low Freq (%.0e Hz)', frequencies(end)), 'Location', 'best');
grid on;
axis equal;

png_file = fullfile(output_dir, 'eis_ref_nyquist.png');
saveas(gcf, png_file);
fprintf('  Nyquist plot: %s\n', png_file);

% Figure 2: Bode plots (reversed frequency axis - high freq on left)
figure(2);
clf;

subplot(2,1,1);
loglog(frequencies, Z_mag_all, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
xlabel('Frequency (Hz)');
ylabel('|Z| (Ohm*m^2)');
title('Bode Plot - Magnitude');
grid on;
set(gca, 'XDir', 'reverse');  % High frequency on left (standard EIS)

subplot(2,1,2);
semilogx(frequencies, phase_deg_all, 'ro-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
title('Bode Plot - Phase');
grid on;
set(gca, 'XDir', 'reverse');  % High frequency on left (standard EIS)

png_file = fullfile(output_dir, 'eis_ref_bode.png');
saveas(gcf, png_file);
fprintf('  Bode plot: %s\n', png_file);

% Figure 3: Least-squares fit quality
figure(3);
clf;
semilogx(frequencies, fit_R2_all, 'ko-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'k');
xlabel('Frequency (Hz)');
ylabel('R^2 (goodness of fit)');
title('Least-Squares Fit Quality');
grid on;
ylim([0.95 1.0]);
set(gca, 'XDir', 'reverse');

png_file = fullfile(output_dir, 'eis_ref_fit_quality.png');
saveas(gcf, png_file);
fprintf('  Fit quality plot: %s\n', png_file);

%% Resistance identification (following paper methodology)
fprintf('\n=== RESISTANCE IDENTIFICATION ===\n');
fprintf('Following Mei et al. (2018) methodology:\n\n');

% Find characteristic points in Nyquist plot
% R_A: Real part at highest frequency (electrode resistance)
R_A = Z_real_all(1);
fprintf('R_A (Point A, electrode resistance):\n');
fprintf('  R_A = %.4e Ohm*m^2 (at %.0e Hz)\n', R_A, frequencies(1));

% R_B: Real part at minimum -Z" (approximately)
[~, idx_min_Zimag] = min(-Z_imag_all);
R_B = Z_real_all(idx_min_Zimag);
fprintf('\nR_B (Point B, end of semicircle):\n');
fprintf('  R_B = %.4e Ohm*m^2 (at %.0e Hz)\n', R_B, frequencies(idx_min_Zimag));

% R_C: Real part at lowest frequency
R_C = Z_real_all(end);
fprintf('\nR_C (Point C, low frequency limit):\n');
fprintf('  R_C = %.4e Ohm*m^2 (at %.0e Hz)\n', R_C, frequencies(end));

% Derived resistances
R_AB = R_B - R_A;  % Semicircle diameter = bulk electrolyte resistance
R_BC = R_C - R_B;  % Warburg/diffusion resistance
R_infinity = R_AB;  % Internal resistance (Re + R∞)

fprintf('\nDerived Resistances:\n');
fprintf('  R_AB = R∞ (bulk electrolyte): %.4e Ohm*m^2\n', R_AB);
fprintf('  R_BC = R_D (diffuse layer): %.4e Ohm*m^2\n', R_BC);
fprintf('  R_∞ (internal resistance): %.4e Ohm*m^2\n', R_infinity);

% Theoretical solution resistance (paper eq. S.21)
F = 96485.3329;  % Faraday constant
R = 8.3144598;   % Gas constant
T = 298.15;      % Temperature (25°C)

% Get actual parameter values from model
DA_val = mphglobal(model, 'DA');
DX_val = mphglobal(model, 'DX');
c_bulk_val = mphglobal(model, 'c_inf');
L_val = mphglobal(model, 'L');
z_val = mphglobal(model, 'z');

% Solution conductivity: κ = (F²/RT) Σ z_i² D_i c_i
% Use element-wise operations (.^, .*, ./) because mphglobal returns arrays
kappa = (F.^2./(R.*T)) .* (z_val.^2 .* DA_val + z_val.^2 .* DX_val) .* c_bulk_val;
R_solution_theory = L_val ./ kappa;

fprintf('\nTheoretical Solution Resistance:\n');
fprintf('  κ (conductivity) = %.4e S/m\n', kappa);
fprintf('  R_s (theory) = L/κ = %.4e Ohm*m^2\n', R_solution_theory);
fprintf('  R_AB (simulation) = %.4e Ohm*m^2\n', R_AB);
fprintf('  Ratio (sim/theory) = %.2f\n', R_AB / R_solution_theory);

%% Display summary
fprintf('\n=== EIS REF SUMMARY ===\n');
fprintf('Implementation: Single electrode blocking electrode (Mei et al. 2018)\n');
fprintf('Impedance method: Least-squares fitting\n');
fprintf('Frequency sweep: %.1e Hz to %.0e Hz (%d points)\n', ...
    frequencies(1), frequencies(end), num_freqs);

fprintf('\nKey Results:\n');
fprintf('  Electrode resistance (R_A): %.4e Ohm*m^2\n', R_A);
fprintf('  Bulk resistance (R_AB): %.4e Ohm*m^2\n', R_AB);
fprintf('  Diffuse layer resistance (R_BC): %.4e Ohm*m^2\n', R_BC);
fprintf('  Fit quality (avg R^2): %.6f\n', mean(fit_R2_all));

fprintf('\n=== DONE ===\n');
fprintf('All results saved in %s/\n', output_dir);

