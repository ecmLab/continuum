%% Debye Length Calculation - COMSOL EDL Benchmark Verification
% This MATLAB script calculates the Debye length using Equation (1) from
% the COMSOL diffuse double layer model documentation
%
% Equation (1): x_D = sqrt(RT*eps_r_D*eps_0 / (2*F^2*c_bulk))

clear; clc;

fprintf('=== DEBYE LENGTH CALCULATION ===\n\n');

%% Physical Constants
R = 8.314;                    % Gas constant [J/(mol·K)]
F = 96485;                    % Faraday constant [C/mol]
eps_0 = 8.854e-12;           % Vacuum permittivity [F/m]

fprintf('Physical Constants:\n');
fprintf('R = %.3f J/(mol·K)\n', R);
fprintf('F = %.0f C/mol\n', F);
fprintf('eps_0 = %.3e F/m\n\n', eps_0);

%% COMSOL Benchmark Parameters
% From diffuse_double_layer_parameters.txt
T0_celsius = 25;             % Temperature [°C]
T0 = T0_celsius + 273.15;    % Temperature [K]
eps_H2O = 78.5;              % Relative permittivity of water
c_bulk = 10;                 % Bulk concentration [mol/m³]
zA = +1;                     % Cation charge
zX = -1;                     % Anion charge

fprintf('COMSOL Benchmark Parameters:\n');
fprintf('T0 = %.1f°C = %.2f K\n', T0_celsius, T0);
fprintf('eps_H2O (eps_r_D) = %.1f\n', eps_H2O);
fprintf('c_bulk = %.0f mol/m³\n', c_bulk);
fprintf('zA = %+d, zX = %+d\n\n', zA, zX);

%% Debye Length Calculation - Equation (1)
% x_D = sqrt(RT*eps_r_D*eps_0 / (2*F^2*c_bulk))

numerator = R * T0 * eps_H2O * eps_0;
denominator = 2 * F^2 * c_bulk;

x_D = sqrt(numerator / denominator);

fprintf('Debye Length Calculation:\n');
fprintf('Numerator = R*T*eps_r_D*eps_0 = %.3e J·F/(mol·m)\n', numerator);
fprintf('Denominator = 2*F²*c_bulk = %.3e C²/(mol²·m³)\n', denominator);
fprintf('x_D = sqrt(%.3e / %.3e)\n', numerator, denominator);
fprintf('x_D = sqrt(%.3e)\n', numerator/denominator);
fprintf('x_D = %.6f m = %.2f nm\n\n', x_D, x_D*1e9);

%% Results Summary
fprintf('=== RESULTS SUMMARY ===\n');
fprintf('Debye Length: x_D = %.2f nm\n', x_D*1e9);

%% Additional EDL Characteristic Lengths
fprintf('\n=== EDL CHARACTERISTIC LENGTHS ===\n');
fprintf('Debye length:        x_D = %.2f nm\n', x_D*1e9);

% Stern layer thickness (from COMSOL parameters)
x_S = 0.2e-9;  % [m]
fprintf('Stern layer:         x_S = %.1f nm\n', x_S*1e9);

% Total EDL thickness (Stern + diffuse layer)
EDL_total = x_S + x_D;
fprintf('Total EDL thickness: %.2f nm\n', EDL_total*1e9);

% COMSOL model domain length
L_cell = x_D * 10;  % From COMSOL parameters
fprintf('Model domain length: L_cell = %.1f nm (10 × x_D)\n', L_cell*1e9);

%% Mesh Size Recommendations (from COMSOL)
h_max = L_cell / 20;        % Maximum element size
h_max_surf = x_D / 100;     % Surface mesh refinement

fprintf('\n=== MESH RECOMMENDATIONS ===\n');
fprintf('Max element size:      h_max = %.2f nm\n', h_max*1e9);
fprintf('Surface refinement:    h_max_surf = %.3f nm\n', h_max_surf*1e9);
fprintf('Surface/bulk ratio:    %.0f:1\n', h_max/h_max_surf);

fprintf('\nCalculation completed successfully!\n');