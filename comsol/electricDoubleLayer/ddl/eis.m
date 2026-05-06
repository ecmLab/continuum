% eis.m
% Electrochemical Impedance Spectroscopy (EIS) - Frequency Sweep
% Based on sinusoidal.m - runs multiple frequencies and collects impedance data
% Created: 2025-10-09

clc; clear;
close all;

%% Create output directory first
output_dir = 'rst/eis';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Setup logging - redirect all output to log file
log_file = fullfile(output_dir, 'eis_simulation.log');
diary(log_file);
diary on;

fprintf('=== EIS SIMULATION LOG ===\n');
fprintf('Started: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('Working directory: %s\n\n', pwd);

%% EIS Configuration
fprintf('=== EIS FREQUENCY SWEEP ===\n\n');

% Frequency range to capture diffusion effects
% With slower diffusion, target 100 Hz to 0.1 Hz range
frequencies = logspace(4, 0.6, 25);  % 100 Hz to 0.1 Hz
num_freqs = length(frequencies);

% Fixed excitation parameters
V_dc = 0.0;               % DC bias voltage (V)
V_ac = 0.010;             % AC amplitude (10 mV)
num_cycles = 10;          % Number of cycles per frequency
pts_per_cycle = 50;       % Time points per cycle

fprintf('EIS Configuration:\n');
fprintf('  Frequencies: ');
for i = 1:num_freqs
    fprintf('%.0e Hz', frequencies(i));
    if i < num_freqs, fprintf(', '); end
end
fprintf('\n');
fprintf('  DC bias: %.3f V\n', V_dc);
fprintf('  AC amplitude: %.3f V (%.1f mV)\n', V_ac, V_ac*1000);
fprintf('  Cycles per frequency: %d\n\n', num_cycles);


%% Initialize storage for EIS results
Z_real_all = zeros(num_freqs, 1);
Z_imag_all = zeros(num_freqs, 1);
Z_mag_all = zeros(num_freqs, 1);
phase_deg_all = zeros(num_freqs, 1);

%% Setup COMSOL (once)
import com.comsol.model.*
import com.comsol.model.util.*

fprintf('Building base COMSOL model...\n');

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

% Set parameters - ASYMMETRIC TWO INTERFACE WITH SINUSOIDAL EXCITATION
model.param.set('T0', '25[degC]', 'Temperature');
model.param.set('V_therm', 'R_const*T0/F_const', 'Thermal voltage');
% Reduce diffusion coefficients to make diffusion process visible in EIS
% (slower diffusion → Warburg impedance appears at measurable frequencies)
model.param.set('DA', '1e-14[m^2/s]', 'Diffusion coefficient, cation');
model.param.set('DX', 'DA*2', 'Diffusion coefficient, anion');
model.param.set('cA_bulk', '10[mol/m^3]', 'Bulk cation concentration');
model.param.set('cX_bulk', 'cA_bulk', 'Bulk anion concentration');
model.param.set('zA', '+1', 'Cation charge');
model.param.set('zX', '-1', 'Anion charge');
model.param.set('Istr_bulk', '0.5*((zA^2+zX^2)*cA_bulk)', 'Bulk ionic strength');
model.param.set('eps_PEO', '10.0', 'Relative permittivity of PEO electrolyte');
model.param.set('xD', 'sqrt(epsilon0_const*eps_PEO*V_therm/(2*F_const*Istr_bulk))', 'Debye length');
model.param.set('xS', '0.2[nm]', 'Stern layer thickness');
model.param.set('L_cell', 'xD*200', 'Total cell length (200 Debye lengths for 2 interfaces)');
model.param.set('h_max', 'L_cell/40', 'Maximum mesh element size');
model.param.set('h_max_surf', 'xD/100', 'Maximum mesh element size (electrodes)');
model.param.set('Cd_GCS', 'epsilon0_const/(xD/eps_PEO+xS/eps_S)', 'Capacitance per unit area (GCS theory)');
model.param.set('eps_S', '10', 'Relative permittivity of Stern layer');

% Sinusoidal voltage parameters (will be updated for each frequency)
model.param.set('V_dc', sprintf('%.6f[V]', V_dc), 'DC bias voltage');
model.param.set('V_ac', sprintf('%.6f[V]', V_ac), 'AC voltage amplitude');
model.param.set('omega', '2*pi*1e5[rad/s]', 'Angular frequency (placeholder)');
model.param.set('phiM_L', 'V_dc + V_ac*sin(omega*t)', 'Left electrode potential (sinusoidal)');
model.param.set('phiM_R', '0[V]', 'Right electrode potential (grounded metal)');

% Create variables for BOTH electrodes
model.component('comp1').variable.create('var1');
% Left electrode variables
model.component('comp1').variable('var1').set('deltaphi_L', 'phiM_L-phi', 'Left electrode-OHP potential difference');
model.component('comp1').variable('var1').set('rho_surf_L', 'epsilon0_const*eps_S*deltaphi_L/xS', 'Left surface charge density');
model.component('comp1').variable('var1').set('i_charging_L', 'd(rho_surf_L,t)', 'Left charging current density (A/m^2)');
% Right electrode variables
model.component('comp1').variable('var1').set('deltaphi_R', 'phiM_R-phi', 'Right electrode-OHP potential difference');
model.component('comp1').variable('var1').set('rho_surf_R', 'epsilon0_const*eps_S*deltaphi_R/xS', 'Right surface charge density');
model.component('comp1').variable('var1').set('i_charging_R', 'd(rho_surf_R,t)', 'Right charging current density (A/m^2)');
% Flux and current variables
model.component('comp1').variable('var1').set('i_cation', 'F_const*zA*tds.tflux_cAx', 'Cation current density (A/m^2)');
model.component('comp1').variable('var1').set('i_anion', 'F_const*zX*tds.tflux_cXx', 'Anion current density (A/m^2)');
model.component('comp1').variable('var1').set('i_total', 'i_cation+i_anion', 'Total ionic current density (A/m^2)');
model.component('comp1').variable('var1').set('V_applied', 'phiM_L', 'Applied voltage at left electrode');

% Create geometry (1D interval from 0 to L_cell)
model.component('comp1').geom('geom1').create('i1', 'Interval');
model.component('comp1').geom('geom1').feature('i1').setIndex('coord', 'L_cell', 1);
model.component('comp1').geom('geom1').run;

% Set physics properties - Electrostatics
model.component('comp1').physics('es').feature('ccn1').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('es').feature('ccn1').set('epsilonr', {'eps_PEO' '0' '0' '0' 'eps_PEO' '0' '0' '0' 'eps_PEO'});

% Boundary conditions - Electrostatics
model.component('comp1').physics('es').create('sfcd1', 'SurfaceChargeDensity', 0);
model.component('comp1').physics('es').feature('sfcd1').selection.set([1]);
model.component('comp1').physics('es').feature('sfcd1').set('rhoqs', 'rho_surf_L');

model.component('comp1').physics('es').create('sfcd2', 'SurfaceChargeDensity', 0);
model.component('comp1').physics('es').feature('sfcd2').selection.set([2]);
model.component('comp1').physics('es').feature('sfcd2').set('rhoqs', 'rho_surf_R');

% Set physics properties - Transport of Diluted Species
model.component('comp1').physics('tds').feature('sp1').setIndex('z', 'zA', 0);
model.component('comp1').physics('tds').feature('sp1').setIndex('z', 'zX', 1);
model.component('comp1').physics('tds').feature('cdm1').set('D_cA', {'DA' '0' '0' '0' 'DA' '0' '0' '0' 'DA'});
model.component('comp1').physics('tds').feature('cdm1').set('D_cX', {'DX' '0' '0' '0' 'DX' '0' '0' '0' 'DX'});
model.component('comp1').physics('tds').feature('init1').setIndex('initc', 'cA_bulk', 0);
model.component('comp1').physics('tds').feature('init1').setIndex('initc', 'cX_bulk', 1);

% Set common properties
model.common('cminpt').set('modified', {'temperature' 'T0'});

% Create mesh with refinement at BOTH electrodes
model.component('comp1').mesh('mesh1').create('edg1', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg1').create('size_L', 'Size');
model.component('comp1').mesh('mesh1').feature('edg1').create('size_R', 'Size');

model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmax', 'h_max');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmaxactive', true);

model.component('comp1').mesh('mesh1').feature('edg1').feature('size_L').selection.geom('geom1', 0);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_L').selection.set([1]);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_L').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_L').set('hmax', 'h_max_surf');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_L').set('hmaxactive', true);

model.component('comp1').mesh('mesh1').feature('edg1').feature('size_R').selection.geom('geom1', 0);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_R').selection.set([2]);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_R').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_R').set('hmax', 'h_max_surf');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size_R').set('hmaxactive', true);

model.component('comp1').mesh('mesh1').run;

fprintf('Base model built successfully!\n');

% Log key model parameters
fprintf('\n=== MODEL PARAMETERS ===\n');
try
    % Get parameter values (some may fail if they're expressions)
    T0_val = model.param.evaluate('T0');
    DA_val = model.param.evaluate('DA');
    DX_val = model.param.evaluate('DX');
    cA_val = model.param.evaluate('cA_bulk');
    xD_val = model.param.evaluate('xD');
    xS_val = model.param.evaluate('xS');
    L_val = model.param.evaluate('L_cell');
    eps_val = model.param.evaluate('eps_PEO');

    fprintf('Temperature: %.2f K\n', T0_val);
    fprintf('Diffusion coefficients:\n');
    fprintf('  D_cation: %.4e m^2/s\n', DA_val);
    fprintf('  D_anion: %.4e m^2/s\n', DX_val);
    fprintf('Bulk concentration: %.4e mol/m^3\n', cA_val);
    fprintf('Debye length: %.4e m\n', xD_val);
    fprintf('Stern layer thickness: %.4e m\n', xS_val);
    fprintf('Domain length: %.4e m\n', L_val);
    fprintf('Permittivity (relative): %.2f\n', eps_val);
catch ME
    fprintf('Warning: Could not evaluate all parameters\n');
    fprintf('  Error: %s\n', ME.message);
end
fprintf('\n');

%% Frequency sweep loop
fprintf('Starting frequency sweep (%d frequencies)...\n', num_freqs);
fprintf('Progress: 0/%d\n\n', num_freqs);

% Will store midpoint info after first simulation
x_mid = [];
mid_idx = [];

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
    fprintf('\n  Running simulation...\n');
    tic;
    model.study('std1').run();
    elapsed = toc;
    fprintf('  Simulation completed in %.2f seconds\n', elapsed);

    %% Get midpoint coordinate (only once, after first simulation)
    if isempty(x_mid)
        x_mid = mphglobal(model, 'L_cell/2');
        x_mid = x_mid(1);

        eval_x = mpheval(model, 'x', 'edim', 1, 'solnum', 1);
        x_coords = eval_x.d1;
        [~, mid_idx] = min(abs(x_coords - x_mid));

        fprintf('  Midpoint: x = %.4e m (mesh index %d)\n', x_mid, mid_idx);
    end

    %% Extract time series data
    fprintf('  Extracting data...\n');

    % Get time and voltage
    time_L = mpheval(model, {'t', 'V_applied'}, 'edim', 0, 'selection', 1);
    t_data = time_L.d1;
    V_applied_t = time_L.d2;

    % Extract i_total at midpoint for all time steps
    num_times = length(t_data);
    i_total_mid_t = zeros(num_times, 1);
    for idx = 1:num_times
        eval_i = mpheval(model, 'i_total', 'edim', 1, 'solnum', idx);
        i_total_mid_t(idx) = eval_i.d1(mid_idx);
    end

    %% Impedance analysis - last 2 cycles
    fprintf('  Analyzing impedance...\n');

    idx_last_cycles = t_data >= (max(t_data) - 2*period);
    t_ss = t_data(idx_last_cycles);
    V_ss = V_applied_t(idx_last_cycles);
    I_ss = i_total_mid_t(idx_last_cycles);

    % Center signals
    V_centered = V_ss - mean(V_ss);
    I_centered = I_ss - mean(I_ss);

    % Apply Hanning window
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
    Fs = 1 / (t_ss(2) - t_ss(1));
    f_vector = Fs * (0:(NFFT/2-1)) / NFFT;

    % Find peak at excitation frequency
    [~, idx_peak] = min(abs(f_vector - freq));

    % Complex impedance
    Z_complex = V_fft(idx_peak) / I_fft(idx_peak);
    Z_mag = abs(Z_complex);
    phase_rad = -angle(Z_complex);
    phase_deg = phase_rad * 180/pi;
    Z_real = real(Z_complex);
    Z_imag = imag(Z_complex);

    % Store results
    Z_real_all(freq_idx) = Z_real;
    Z_imag_all(freq_idx) = Z_imag;
    Z_mag_all(freq_idx) = Z_mag;
    phase_deg_all(freq_idx) = phase_deg;

    fprintf('\n  Results:\n');
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

% Save CSV
csv_file = fullfile(output_dir, 'eis_data.csv');
eis_table = table(frequencies', Z_real_all, Z_imag_all, Z_mag_all, phase_deg_all, ...
    'VariableNames', {'freq_Hz', 'Z_real_Ohm_m2', 'Z_imag_Ohm_m2', 'Z_mag_Ohm_m2', 'phase_deg'});
writetable(eis_table, csv_file);
fprintf('  CSV: %s\n', csv_file);

% Save MAT file
mat_file = fullfile(output_dir, 'eis_results.mat');
save(mat_file, 'frequencies', 'Z_real_all', 'Z_imag_all', 'Z_mag_all', 'phase_deg_all', ...
    'V_dc', 'V_ac', 'num_cycles');
fprintf('  MAT: %s\n', mat_file);

% Save COMSOL model with solution
fprintf('  Saving COMSOL model...\n');
mph_file = fullfile(output_dir, 'eis_solved.mph');
mphsave(model, mph_file);
fprintf('  COMSOL model: %s\n', mph_file);

%% Generate plots
fprintf('\nGenerating plots...\n');

% Figure 1: Nyquist plot
figure(1);
plot(Z_real_all, -Z_imag_all, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;
plot(Z_real_all(1), -Z_imag_all(1), 'rs', 'MarkerSize', 12, 'LineWidth', 2);
plot(Z_real_all(end), -Z_imag_all(end), 'gs', 'MarkerSize', 12, 'LineWidth', 2);
xlabel('Z'' (Ohm*m^2)');
ylabel('-Z" (Ohm*m^2)');
title('Nyquist Plot (EIS)');
legend('Data', sprintf('High (%.0e Hz)', frequencies(1)), ...
    sprintf('Low (%.0e Hz)', frequencies(end)), 'Location', 'best');
grid on;
xlim([0.15 0.2]);
saveas(gcf, fullfile(output_dir, 'eis_nyquist.png'));

% Figure 2: Bode plots
figure(2);
subplot(2,1,1);
loglog(frequencies, Z_mag_all, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
xlabel('Frequency (Hz)');
ylabel('|Z| (Ohm*m^2)');
title('Bode Plot - Magnitude');
grid on;
set(gca, 'XDir', 'reverse');  % High frequency on left

subplot(2,1,2);
semilogx(frequencies, phase_deg_all, 'ro-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
title('Bode Plot - Phase');
grid on;
ylim([-5 100]);
set(gca, 'XDir', 'reverse');  % High frequency on left

saveas(gcf, fullfile(output_dir, 'eis_bode.png'));

fprintf('  Plots saved to %s/\n', output_dir);

%% Display summary
fprintf('\n=== EIS SUMMARY ===\n');
fprintf('Frequency sweep: %.0e Hz to %.0e Hz (%d points)\n', ...
    frequencies(1), frequencies(end), num_freqs);
fprintf('\nImpedance results:\n');
for i = 1:num_freqs
    fprintf('  %.0e Hz: Z = %.4e Ohm*m^2, phi = %.2f deg\n', ...
        frequencies(i), Z_mag_all(i), phase_deg_all(i));
end

fprintf('\n=== DONE ===\n');
fprintf('Completed: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('All results saved in %s/\n', output_dir);
fprintf('\nOutput files:\n');
fprintf('  - %s (CSV data)\n', csv_file);
fprintf('  - %s (MATLAB data)\n', mat_file);
fprintf('  - %s (COMSOL model)\n', mph_file);
fprintf('  - %s (simulation log)\n', log_file);
fprintf('  - %s/eis_nyquist.png (Nyquist plot)\n', output_dir);
fprintf('  - %s/eis_bode.png (Bode plot)\n', output_dir);

% Close diary
diary off;
fprintf('\nLog file saved: %s\n', log_file);
