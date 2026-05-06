% oneInterface_v1.m
% Script to run the diffuse double layer COMSOL model for single electrode
% Created: 2025-10-08
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

fprintf('Building model from scratch...\n');

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

% Set parameters (aligned with Mei et al. 2018 JPCC paper)
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
model.param.set('L_cell', '160[nm]', 'Cell length (electrolyte domain thickness from paper Table 2)');
model.param.set('Le', '10[nm]', 'Electrode thickness (from paper Table 1)');

% Mesh parameters
model.param.set('h_max', 'L_cell/40', 'Maximum mesh element size');
model.param.set('h_max_surf', 'xD/100', 'Maximum mesh element size (electrode)');

% Capacitance (for reference)
model.param.set('Cd_GCS', 'epsilon0_const/(xD/eps_r+xS/eps_r)', 'Capacitance per unit area (GCS theory)');

% Applied potential
model.param.set('phiM', '0.01[V]', 'DC potential at electrode (psi_dc from paper)');

% Create variables
model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('deltaphi', 'phiM-phi', 'Electrode-OHP potential difference');
model.component('comp1').variable('var1').set('rho_surf', 'epsilon0_const*eps_r*deltaphi/xS', 'Surface charge density');
model.component('comp1').variable('var1').set('i_charging', 'd(rho_surf,t)', 'Charging current density from surface charge (A/m^2)');
model.component('comp1').variable('var1').set('i_cation', 'F_const*zA*tds.tflux_cAx', 'Cation current density (A/m^2)');
model.component('comp1').variable('var1').set('i_anion', 'F_const*zX*tds.tflux_cXx', 'Anion current density (A/m^2)');
model.component('comp1').variable('var1').set('i_total', 'i_cation+i_anion', 'Total ionic current density (A/m^2)');
% Try alternative flux definitions
model.component('comp1').variable('var1').set('flux_cA_diffusive', 'tds.dflux_cAx', 'Cation diffusive flux (mol/m^2/s)');
model.component('comp1').variable('var1').set('flux_cA_migration', 'tds.mflux_cAx', 'Cation migration flux (mol/m^2/s)');
model.component('comp1').variable('var1').set('flux_cX_diffusive', 'tds.dflux_cXx', 'Anion diffusive flux (mol/m^2/s)');
model.component('comp1').variable('var1').set('flux_cX_migration', 'tds.mflux_cXx', 'Anion migration flux (mol/m^2/s)');
% Manual flux calculation from gradients
model.component('comp1').variable('var1').set('flux_cA_diff_manual', '-DA*cAx', 'Cation diffusive flux (manual, mol/m^2/s)');
model.component('comp1').variable('var1').set('flux_cA_mig_manual', '-DA*zA*cA*F_const/(R_const*T0)*phix', 'Cation migration flux (manual, mol/m^2/s)');
model.component('comp1').variable('var1').set('E_field', '-phix', 'Electric field (V/m)');

% Create geometry (1D interval)
model.component('comp1').geom('geom1').create('i1', 'Interval');
model.component('comp1').geom('geom1').feature('i1').setIndex('coord', 'L_cell', 1);
model.component('comp1').geom('geom1').run;

% Set physics properties - Electrostatics
model.component('comp1').physics('es').feature('ccn1').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('es').feature('ccn1').set('epsilonr', {'eps_r' '0' '0' '0' 'eps_r' '0' '0' '0' 'eps_r'});
model.component('comp1').physics('es').create('sfcd1', 'SurfaceChargeDensity', 0);
model.component('comp1').physics('es').feature('sfcd1').selection.set([1]);
model.component('comp1').physics('es').feature('sfcd1').set('rhoqs', 'rho_surf');
model.component('comp1').physics('es').create('gnd1', 'Ground', 0);
model.component('comp1').physics('es').feature('gnd1').selection.set([2]);

% Set physics properties - Transport of Diluted Species
model.component('comp1').physics('tds').feature('sp1').setIndex('z', 'zA', 0);
model.component('comp1').physics('tds').feature('sp1').setIndex('z', 'zX', 1);
model.component('comp1').physics('tds').feature('cdm1').set('D_cA', {'DA' '0' '0' '0' 'DA' '0' '0' '0' 'DA'});
model.component('comp1').physics('tds').feature('cdm1').set('D_cX', {'DX' '0' '0' '0' 'DX' '0' '0' '0' 'DX'});
model.component('comp1').physics('tds').feature('init1').setIndex('initc', 'cA_bulk', 0);
model.component('comp1').physics('tds').feature('init1').setIndex('initc', 'cX_bulk', 1);
model.component('comp1').physics('tds').create('conc1', 'Concentration', 0);
model.component('comp1').physics('tds').feature('conc1').selection.set([2]);
model.component('comp1').physics('tds').feature('conc1').setIndex('species', true, 0);
model.component('comp1').physics('tds').feature('conc1').setIndex('c0', 'cA_bulk', 0);
model.component('comp1').physics('tds').feature('conc1').setIndex('species', true, 1);
model.component('comp1').physics('tds').feature('conc1').setIndex('c0', 'cX_bulk', 1);

% Set common properties
model.common('cminpt').set('modified', {'temperature' 'T0'});

% Create mesh with refined electrode region (matching ddl0.m)
model.component('comp1').mesh('mesh1').create('edg1', 'Edge');
model.component('comp1').mesh('mesh1').feature('edg1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('edg1').create('size2', 'Size');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size2').selection.geom('geom1', 0);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size2').selection.set([1]);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmax', 'h_max');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('edg1').feature('size2').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size2').set('hmax', 'h_max_surf');
model.component('comp1').mesh('mesh1').feature('edg1').feature('size2').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').run;

% Configure study settings (matching ddl0.m)
% Stationary study is inactive (only for parametric sweep if needed)
model.study('std1').feature('stat').active(false);

% Time-dependent study settings
model.study('std1').feature('time').set('tunit', 'ms');
% Use finer steps under 1 ms while keeping the original span to 1000 ms (double resolution <0.01 ms)
model.study('std1').feature('time').set('tlist', '10^{range(log10(0.001),1/12,log10(0.01))} 10^{range(log10(0.01)+1/6,1/6,log10(1))} 10^{range(log10(1),1/3,log10(1000))}');

fprintf('Model built successfully!\n');
fprintf('  Model tag: %s\n', char(model.tag()));
fprintf('  Components: %d\n', model.component.size());
fprintf('  Physics modules: %d\n', model.component('comp1').physics.size());

%% Run the time-dependent study
fprintf('\nRunning time-dependent simulation...\n');
fprintf('Time range: 0.001 to 1000 ms (logarithmic)\n');
tic;
model.study('std1').run();
elapsed_time = toc;
fprintf('Simulation completed in %.2f seconds\n', elapsed_time);

%% Extract results along the domain (final time step)
fprintf('\nExtracting results at final time step...\n');

% Use mpheval to extract data along the geometry at last time step
% For 1D geometry, extract along the edge
eval_result = mpheval(model, {'x', 'phi', 'cA', 'cX', 'es.normE', 'es.rhoq'}, 'edim', 1, 'solnum', 'end');

% Extract data from result structure
x_data = eval_result.d1;         % x coordinates
phi_data = eval_result.d2;       % electric potential
cA_data = eval_result.d3;        % cation concentration
cX_data = eval_result.d4;        % anion concentration
Es_data = eval_result.d5;        % electric field magnitude
rho_data = eval_result.d6;       % space charge density

fprintf('Data extracted successfully!\n');
fprintf('  Number of spatial points: %d\n', length(x_data));

% Extract time series data of potential and concentrations at electrode position (left boundary, x=0)
fprintf('\nExtracting time evolution at electrode position (x=0)...\n');
time_result = mpheval(model, {'t', 'phi', 'cA', 'cX'}, 'edim', 0, 'selection', 1);
t_data = time_result.d1;         % time points
phi_t = time_result.d2;          % potential vs time at bulk
cA_t = time_result.d3;           % cation concentration vs time
cX_t = time_result.d4;           % anion concentration vs time

% Extract time series data of current density at bulk position (right boundary, x=L_cell)
fprintf('\nExtracting time evolution at bulk position (x=L_cell)...\n');
time_crnt = mpheval(model, {'i_cation', 'i_anion', 'i_total'}, 'edim', 0, 'selection', 2);
i_cation_t = time_crnt.d1;     % cation current density vs time
i_anion_t = time_crnt.d2;      % anion current density vs time
i_total_t = time_crnt.d3;      % total current density vs time

fprintf('  Number of time points: %d\n', length(t_data));

%% Create output directory
output_dir = 'rst/oneInterface';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Created output directory: %s\n', output_dir);
end

%% Add Results plots to COMSOL model for visualization (matching PNG output figures)
fprintf('\nAdding result plots to COMSOL model...\n');

% Figure 1: Electric Potential Distribution (Final State) - oneInterface_potential.png
model.result.create('pg1', 'PlotGroup1D');
model.result('pg1').set('data', 'dset1');
model.result('pg1').label('Figure 1: Electric Potential Distribution (Final State)');
model.result('pg1').set('xlabel', 'Distance from electrode (nm)');
model.result('pg1').set('ylabel', 'Electric potential (mV)');
model.result('pg1').set('xlabelactive', true);
model.result('pg1').set('ylabelactive', true);
model.result('pg1').set('titletype', 'label');
model.result('pg1').create('lngr1', 'LineGraph');
model.result('pg1').feature('lngr1').selection.set([1]);
model.result('pg1').feature('lngr1').set('expr', 'phi');
model.result('pg1').feature('lngr1').set('unit', 'mV');
model.result('pg1').feature('lngr1').set('xdata', 'expr');
model.result('pg1').feature('lngr1').set('xdataexpr', 'x');
model.result('pg1').feature('lngr1').set('xdataunit', 'nm');
model.result('pg1').feature('lngr1').set('resolution', 'normal');

% Figure 2: Ion Concentration Profiles (Final State) - oneInterface_concentrations.png
model.result.create('pg2', 'PlotGroup1D');
model.result('pg2').set('data', 'dset1');
model.result('pg2').label('Figure 2: Ion Concentration Profiles (Final State)');
model.result('pg2').set('xlabel', 'Distance from electrode (nm)');
model.result('pg2').set('ylabel', 'Concentration (mol/m^3)');
model.result('pg2').set('xlabelactive', true);
model.result('pg2').set('ylabelactive', true);
model.result('pg2').set('titletype', 'label');
model.result('pg2').create('lngr1', 'LineGraph');
model.result('pg2').feature('lngr1').selection.set([1]);
model.result('pg2').feature('lngr1').set('expr', 'cA');
model.result('pg2').feature('lngr1').set('unit', 'mol/m^3');
model.result('pg2').feature('lngr1').set('descr', 'Cation (A^+)');
model.result('pg2').feature('lngr1').set('xdata', 'expr');
model.result('pg2').feature('lngr1').set('xdataexpr', 'x');
model.result('pg2').feature('lngr1').set('xdataunit', 'nm');
model.result('pg2').feature('lngr1').set('resolution', 'normal');
model.result('pg2').create('lngr2', 'LineGraph');
model.result('pg2').feature('lngr2').selection.set([1]);
model.result('pg2').feature('lngr2').set('expr', 'cX');
model.result('pg2').feature('lngr2').set('unit', 'mol/m^3');
model.result('pg2').feature('lngr2').set('descr', 'Anion (X^-)');
model.result('pg2').feature('lngr2').set('xdata', 'expr');
model.result('pg2').feature('lngr2').set('xdataexpr', 'x');
model.result('pg2').feature('lngr2').set('xdataunit', 'nm');
model.result('pg2').feature('lngr2').set('resolution', 'normal');

% Figure 3: Electric Field Distribution (Final State) - oneInterface_efield.png
model.result.create('pg3', 'PlotGroup1D');
model.result('pg3').set('data', 'dset1');
model.result('pg3').label('Figure 3: Electric Field Distribution (Final State)');
model.result('pg3').set('xlabel', 'Distance from electrode (nm)');
model.result('pg3').set('ylabel', 'Electric field (MV/m)');
model.result('pg3').set('xlabelactive', true);
model.result('pg3').set('ylabelactive', true);
model.result('pg3').set('titletype', 'label');
model.result('pg3').create('lngr1', 'LineGraph');
model.result('pg3').feature('lngr1').selection.set([1]);
model.result('pg3').feature('lngr1').set('expr', 'E_field');
model.result('pg3').feature('lngr1').set('unit', 'MV/m');
model.result('pg3').feature('lngr1').set('xdata', 'expr');
model.result('pg3').feature('lngr1').set('xdataexpr', 'x');
model.result('pg3').feature('lngr1').set('xdataunit', 'nm');
model.result('pg3').feature('lngr1').set('resolution', 'normal');

% Figure 4: Space Charge Density (Final State) - oneInterface_charge.png
model.result.create('pg4', 'PlotGroup1D');
model.result('pg4').set('data', 'dset1');
model.result('pg4').label('Figure 4: Space Charge Density (Final State)');
model.result('pg4').set('xlabel', 'Distance from electrode (nm)');
model.result('pg4').set('ylabel', 'Charge density (C/m^3)');
model.result('pg4').set('xlabelactive', true);
model.result('pg4').set('ylabelactive', true);
model.result('pg4').set('titletype', 'label');
model.result('pg4').create('lngr1', 'LineGraph');
model.result('pg4').feature('lngr1').selection.set([1]);
model.result('pg4').feature('lngr1').set('expr', 'es.rhoq');
model.result('pg4').feature('lngr1').set('unit', 'C/m^3');
model.result('pg4').feature('lngr1').set('xdata', 'expr');
model.result('pg4').feature('lngr1').set('xdataexpr', 'x');
model.result('pg4').feature('lngr1').set('xdataunit', 'nm');
model.result('pg4').feature('lngr1').set('resolution', 'normal');

% Figure 5: Time Evolution at Electrode - oneInterface_time_evolution.png
model.result.create('pg5', 'PlotGroup1D');
model.result('pg5').set('data', 'dset1');
model.result('pg5').label('Figure 5: Time Evolution at Electrode');
model.result('pg5').set('xlabel', 'Time (ms)');
model.result('pg5').set('ylabel', 'Values');
model.result('pg5').set('xlabelactive', true);
model.result('pg5').set('ylabelactive', false);
model.result('pg5').set('titletype', 'label');
model.result('pg5').set('axislimits', 'on');
model.result('pg5').set('xlog', true);
% Subplot 1: Potential
model.result('pg5').create('ptgr1', 'PointGraph');
model.result('pg5').feature('ptgr1').selection.set([1]);  % Left electrode (boundary 1)
model.result('pg5').feature('ptgr1').set('expr', 'phi');
model.result('pg5').feature('ptgr1').set('unit', 'mV');
model.result('pg5').feature('ptgr1').set('descr', 'Potential at electrode');
model.result('pg5').feature('ptgr1').set('xdata', 'expr');
model.result('pg5').feature('ptgr1').set('xdataexpr', 't');
model.result('pg5').feature('ptgr1').set('xdataunit', 'ms');
% Subplot 2: Cation concentration
model.result('pg5').create('ptgr2', 'PointGraph');
model.result('pg5').feature('ptgr2').selection.set([1]);
model.result('pg5').feature('ptgr2').set('expr', 'cA');
model.result('pg5').feature('ptgr2').set('unit', 'mol/m^3');
model.result('pg5').feature('ptgr2').set('descr', 'Cation concentration at electrode');
model.result('pg5').feature('ptgr2').set('xdata', 'expr');
model.result('pg5').feature('ptgr2').set('xdataexpr', 't');
model.result('pg5').feature('ptgr2').set('xdataunit', 'ms');
% Subplot 3: Anion concentration
model.result('pg5').create('ptgr3', 'PointGraph');
model.result('pg5').feature('ptgr3').selection.set([1]);
model.result('pg5').feature('ptgr3').set('expr', 'cX');
model.result('pg5').feature('ptgr3').set('unit', 'mol/m^3');
model.result('pg5').feature('ptgr3').set('descr', 'Anion concentration at electrode');
model.result('pg5').feature('ptgr3').set('xdata', 'expr');
model.result('pg5').feature('ptgr3').set('xdataexpr', 't');
model.result('pg5').feature('ptgr3').set('xdataunit', 'ms');

% Figure 6: Current Density Evolution - oneInterface_current_evolution.png
model.result.create('pg6', 'PlotGroup1D');
model.result('pg6').set('data', 'dset1');
model.result('pg6').label('Figure 6: Current Density Evolution During Double Layer Charging');
model.result('pg6').set('xlabel', 'Time (ms)');
model.result('pg6').set('ylabel', 'Current Density (A/m^2)');
model.result('pg6').set('xlabelactive', true);
model.result('pg6').set('ylabelactive', true);
model.result('pg6').set('titletype', 'label');
model.result('pg6').set('xlog', true);
% Total current
model.result('pg6').create('ptgr1', 'PointGraph');
model.result('pg6').feature('ptgr1').selection.set([2]);  % Right bulk (boundary 2)
model.result('pg6').feature('ptgr1').set('expr', 'i_total');
model.result('pg6').feature('ptgr1').set('unit', 'A/m^2');
model.result('pg6').feature('ptgr1').set('descr', 'Total Current');
model.result('pg6').feature('ptgr1').set('xdata', 'expr');
model.result('pg6').feature('ptgr1').set('xdataexpr', 't');
model.result('pg6').feature('ptgr1').set('xdataunit', 'ms');
% Cation contribution
model.result('pg6').create('ptgr2', 'PointGraph');
model.result('pg6').feature('ptgr2').selection.set([2]);
model.result('pg6').feature('ptgr2').set('expr', 'i_cation');
model.result('pg6').feature('ptgr2').set('unit', 'A/m^2');
model.result('pg6').feature('ptgr2').set('descr', 'Cation Contribution');
model.result('pg6').feature('ptgr2').set('xdata', 'expr');
model.result('pg6').feature('ptgr2').set('xdataexpr', 't');
model.result('pg6').feature('ptgr2').set('xdataunit', 'ms');
% Anion contribution
model.result('pg6').create('ptgr3', 'PointGraph');
model.result('pg6').feature('ptgr3').selection.set([2]);
model.result('pg6').feature('ptgr3').set('expr', 'i_anion');
model.result('pg6').feature('ptgr3').set('unit', 'A/m^2');
model.result('pg6').feature('ptgr3').set('descr', 'Anion Contribution');
model.result('pg6').feature('ptgr3').set('xdata', 'expr');
model.result('pg6').feature('ptgr3').set('xdataexpr', 't');
model.result('pg6').feature('ptgr3').set('xdataunit', 'ms');

fprintf('Result plots added to model (matching PNG output figures):\n');
fprintf('  - pg1: Figure 1 - Electric Potential Distribution (Final State)\n');
fprintf('  - pg2: Figure 2 - Ion Concentration Profiles (Final State)\n');
fprintf('  - pg3: Figure 3 - Electric Field Distribution (Final State)\n');
fprintf('  - pg4: Figure 4 - Space Charge Density (Final State)\n');
fprintf('  - pg5: Figure 5 - Time Evolution at Electrode\n');
fprintf('  - pg6: Figure 6 - Current Density Evolution\n');

%% Remove old plot groups if they exist (pg7, pg8)
try
    model.result.remove('pg7');
    model.result.remove('pg8');
catch
    % Plots don't exist, continue
end

%% Save COMSOL model file
model_file = fullfile(output_dir, 'oneInterface.mph');
mphsave(model, model_file);
fprintf('\nCOMSOL model saved to: %s\n', model_file);
fprintf('All 6 plot groups (pg1-pg6) are included in the .mph file\n');

%% Save results to file
results_file = fullfile(output_dir, 'oneInterface_results.mat');
save(results_file, 'x_data', 'phi_data', 'cA_data', 'cX_data', 'Es_data', 'rho_data', ...
    't_data', 'phi_t', 'cA_t', 'cX_t', 'i_cation_t', 'i_anion_t', 'i_total_t', ...
    'elapsed_time');
fprintf('Results saved to: %s\n', results_file);

%% Save spatial data as CSV (final time step)
csv_file = fullfile(output_dir, 'oneInterface_data_spatial.csv');
data_table = table(x_data*1e9, phi_data*1000, cA_data, cX_data, Es_data/1e6, rho_data, ...
    'VariableNames', {'x_nm', 'phi_mV', 'cA_molm3', 'cX_molm3', 'E_MVm', 'rho_Cm3'});
writetable(data_table, csv_file);
fprintf('Spatial data saved to CSV: %s\n', csv_file);

%% Save time evolution data as CSV (bulk position)
csv_time_file = fullfile(output_dir, 'oneInterface_data_time.csv');
time_table = table(t_data*1000, phi_t*1000, cA_t, cX_t, i_cation_t, i_anion_t, i_total_t, ...
    'VariableNames', {'t_ms', 'phi_bulk_mV', 'cA_bulk_molm3', 'cX_bulk_molm3', ...
    'i_cation_Am2', 'i_anion_Am2', 'i_total_Am2'});
writetable(time_table, csv_time_file);
fprintf('Time evolution data (bulk) saved to CSV: %s\n', csv_time_file);

%% Create plots
fprintf('\nGenerating plots...\n');

% Figure 1: Electric potential (final time)
figure(1);
plot(x_data*1e9, phi_data*1000, 'b-', 'LineWidth', 2);
xlabel('Distance from electrode (nm)');
ylabel('Electric potential (mV)');
title('Electric Potential Distribution (Final State)');
grid on;
saveas(gcf, fullfile(output_dir, 'oneInterface_potential.png'));

% Figure 2: Ion concentrations (final time)
figure(2);
plot(x_data*1e9, cA_data, 'r-', 'LineWidth', 2); hold on;
plot(x_data*1e9, cX_data, 'b-', 'LineWidth', 2);
xlabel('Distance from electrode (nm)');
ylabel('Concentration (mol/m^3)');
legend('Cation (A^+)', 'Anion (X^-)');
title('Ion Concentration Profiles (Final State)');
grid on;
saveas(gcf, fullfile(output_dir, 'oneInterface_concentrations.png'));

% Figure 3: Electric field (final time)
figure(3);
plot(x_data*1e9, Es_data/1e6, 'g-', 'LineWidth', 2);
xlabel('Distance from electrode (nm)');
ylabel('Electric field (MV/m)');
title('Electric Field Distribution (Final State)');
grid on;
saveas(gcf, fullfile(output_dir, 'oneInterface_efield.png'));

% Figure 4: Space charge density (final time)
figure(4);
plot(x_data*1e9, rho_data, 'm-', 'LineWidth', 2);
xlabel('Distance from electrode (nm)');
ylabel('Charge density (C/m^3)');
title('Space Charge Density (Final State)');
grid on;
saveas(gcf, fullfile(output_dir, 'oneInterface_charge.png'));

% Figure 5: Time evolution at electrode
figure(5);
subplot(3,1,1);
semilogx(t_data*1000, phi_t*1000, 'b-', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Potential (mV)');
title('Electric Potential at Electrode vs Time');
grid on;

subplot(3,1,2);
semilogx(t_data*1000, cA_t, 'r-', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Concentration (mol/m^3)');
title('Cation Concentration at Electrode vs Time');
grid on;

subplot(3,1,3);
semilogx(t_data*1000, cX_t, 'b-', 'LineWidth', 2);
xlabel('Time (ms)');
ylabel('Concentration (mol/m^3)');
title('Anion Concentration at Electrode vs Time');
grid on;

saveas(gcf, fullfile(output_dir, 'oneInterface_time_evolution.png'));

% Figure 6: Current density vs time (separate detailed plot)
figure(6);
semilogx(t_data*1000, i_total_t, 'k-', 'LineWidth', 2); hold on;
semilogx(t_data*1000, i_cation_t, 'r--', 'LineWidth', 1.5);
semilogx(t_data*1000, i_anion_t, 'b--', 'LineWidth', 1.5);
xlabel('Time (ms)');
ylabel('Current Density (A/m^2)');
title('Current Density Evolution During Double Layer Charging');
legend('Total Current', 'Cation Contribution', 'Anion Contribution', 'Location', 'best');
grid on;
saveas(gcf, fullfile(output_dir, 'oneInterface_current_evolution.png'));

fprintf('Plots saved to %s/\n', output_dir);

%% Display key parameters
fprintf('\n=== Model Parameters ===\n');
fprintf('Temperature: %s\n', char(model.param.get('T0')));
fprintf('Thermal voltage: %s\n', char(model.param.get('V_therm')));
fprintf('Diffusion coefficient: %s\n', char(model.param.get('DA')));
fprintf('Bulk concentration: %s\n', char(model.param.get('cA_bulk')));
fprintf('Debye length: %s\n', char(model.param.get('xD')));
fprintf('Stern layer thickness: %s\n', char(model.param.get('xS')));
fprintf('Electrode potential vs PZC: %s\n', char(model.param.get('phiM')));

%% Display key results
fprintf('\n=== Key Results ===\n');
fprintf('Electrode potential: %.4f mV\n', phi_data(1)*1000);
fprintf('Bulk potential: %.4f mV\n', phi_data(end)*1000);
fprintf('Potential drop: %.4f mV\n', (phi_data(1)-phi_data(end))*1000);
fprintf('Max cation concentration: %.2f mol/m^3\n', max(cA_data));
fprintf('Min anion concentration: %.2f mol/m^3\n', min(cX_data));
fprintf('Max charge density: %.2e C/m^3\n', max(abs(rho_data)));
fprintf('Max electric field: %.2f MV/m\n', max(Es_data)/1e6);

fprintf('\n=== DONE ===\n');
fprintf('All results saved in ./%s/ directory:\n', output_dir);
fprintf('  - ddl_solved.mph (COMSOL model with time-dependent solution)\n');
fprintf('  - ddl_results.mat (MATLAB data with spatial and temporal data)\n');
fprintf('  - ddl_data_spatial.csv (Spatial profiles at final time)\n');
fprintf('  - ddl_data_time.csv (Time evolution: phi, cA, cX, currents)\n');
fprintf('  - ddl_potential.png (Spatial distribution)\n');
fprintf('  - ddl_concentrations.png (Spatial distribution)\n');
fprintf('  - ddl_efield.png (Spatial distribution)\n');
fprintf('  - ddl_charge.png (Spatial distribution)\n');
fprintf('  - ddl_time_evolution.png (Time evolution: phi, cA, cX, currents)\n');
fprintf('  - ddl_current_evolution.png (Current density vs time)\n');
