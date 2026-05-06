%{ 
Dan Brunner
CEWLab
Stability Diagram (Damping)
%}
tic

clear
clc

%% Display Current Time
currentTime = datetime('now');
disp(['Current Time: ', datestr(currentTime)]);

%% Inputs

% % Li+
% V_AC = 1.4339e+11 * 2 / 2; % [Vpeak]; 2V_AC = AC Voltage [Vpeak]; 4V_AC = AC Voltage [Vpeak-to-peak]
% U_DC = 41.2251 * 2 / 2; % [V]
% % r_0 = 5.531e-3; % [m]
% r_0 = 5e-3; % [m]
% f = 2e5; % [Hz] NOTE: decreasing magnitude by one AND increasing magnitude by one both result in instability; there is an optimal f for a given k
% atomic_numbers = [3];
% [q,a,k,omega,Q_zeta_ratio] = Parameters(V_AC,U_DC,r_0,f,atomic_numbers);
% q = 2e11;
% a = 140; % Stable around 110-115 for r=30000

% % Na+
% V_AC = 1.4339e+11 * 2 / 2; % [Vpeak]; 2V_AC = AC Voltage [Vpeak]; 4V_AC = AC Voltage [Vpeak-to-peak]
% U_DC = 41.2251 * 2 / 2; % [V]
% % r_0 = 5.531e-3; % [m]
% r_0 = 5e-3; % [m]
% f = 2e5; % [Hz] NOTE: decreasing magnitude by one AND increasing magnitude by one both result in instability; there is an optimal f for a given k
% atomic_numbers = [11];
% [q,a,k,omega,Q_zeta_ratio] = Parameters(V_AC,U_DC,r_0,f,atomic_numbers);
% q = 2e11 * 1/3.284431784798363;
% a = 140 * 1/3.284431784798363; 

% % 2micron microsphere
% % V_AC = 40 / 2; % [Vpeak]; 2V_AC = AC Voltage [Vpeak]; 4V_AC = AC Voltage [Vpeak-to-peak]
% % U_DC = 0.002 / 2; % [V]
% V_AC = 10 / 2; % [Vpeak]; 2V_AC = AC Voltage [Vpeak]; 4V_AC = AC Voltage [Vpeak-to-peak]
% U_DC = 0.002 / 2; % [V]
% % V_AC = 5 / 2; % [Vpeak]; 2V_AC = AC Voltage [Vpeak]; 4V_AC = AC Voltage [Vpeak-to-peak]
% % U_DC = 0.002 / 2; % [V]
% r_0 = 5.531e-3; % [m]
% f = 1e4; % [Hz] NOTE: decreasing magnitude by one AND increasing magnitude by one both result in instability; there is an optimal f for a given k
% atomic_numbers = [128]; % 124-128
% [q,a,k,omega,Q_zeta_ratio] = Parameters(V_AC,U_DC,r_0,f,atomic_numbers);

% 10micron microsphere
% V_AC = 40 / 2; % [Vpeak]; 2V_AC = AC Voltage [Vpeak]; 4V_AC = AC Voltage [Vpeak-to-peak]
% U_DC = 0.002 / 2; % [V]
V_AC = 10 / 2; % [Vpeak]; 2V_AC = AC Voltage [Vpeak]; 4V_AC = AC Voltage [Vpeak-to-peak]
U_DC = 0.002 / 2; % [V]
% V_AC = 5 / 2; % [Vpeak]; 2V_AC = AC Voltage [Vpeak]; 4V_AC = AC Voltage [Vpeak-to-peak]
% U_DC = 0.002 / 2; % [V]
r_0 = 5.531e-3; % [m]
f = 1e4; % [Hz] NOTE: decreasing magnitude by one AND increasing magnitude by one both result in instability; there is an optimal f for a given k
atomic_numbers = [119]; % 119-123
[q,a,k,omega,Q_zeta_ratio] = Parameters(V_AC,U_DC,r_0,f,atomic_numbers);

%% Point Calculation

r_max = 30000; % 30000 (approximate max) takes ~915s (15min 15s)
r_min = -r_max;
r_points = r_max - r_min + 1;
r = linspace(r_min,r_max,r_points);
               
% Initialize the infinite determinant matrix
inf_matrix_zero_x = eye(r_points);
inf_matrix_zero_y = eye(r_points);

% Calculate all of the xi's necessary for the infinite determinant matrix
xi_x = q ./ (4 * r.^2 - a + k^2 / 4);
xi_y = q ./ (4 * r.^2 + a + k^2 / 4);
        
for f = 1:r_points - 1
    y(f) = f + 1;
    x(f) = f - 1;
    
    % Left diagonal of infinite matrix                
    inf_matrix_zero_x(y(f),f) = xi_x(y(f));
    inf_matrix_zero_y(y(f),f) = xi_y(y(f));
    
    % Right diagonal of infinite matrix
    inf_matrix_zero_x(f,y(f)) = xi_x(f);
    inf_matrix_zero_y(f,y(f)) = xi_y(f);
end
        
% Calculate the determinant of the infinite matrix and its derivative
inf_det_zero_x = det(inf_matrix_zero_x);
inf_det_zero_y = det(inf_matrix_zero_y);

% Calculate alpha

if single(k) > 380
    alpha_x = 4 / pi * asinh(sqrt(-inf_det_zero_x * sin(sym(pi / 2 * sqrt(a - k^2 / 4)))^2));
    alpha_y = 4 / pi * asinh(sqrt(-inf_det_zero_y * sin(sym(pi / 2 * sqrt(-a - k^2 / 4)))^2));
else
    alpha_x = 4 / pi * asinh(sqrt(-inf_det_zero_x .* sin(pi / 2 * sqrt(a - k^2 / 4))^2));
    alpha_y = 4 / pi * asinh(sqrt(-inf_det_zero_y .* sin(pi / 2 * sqrt(-a - k^2 / 4))^2));
end

% clearvars -except mu k

if single(alpha_x < k) && single(alpha_y < k)
    disp('Stable')
elseif single(alpha_x > k) && single(alpha_y > k)
    disp('Unstable')
elseif single(alpha_x < k) && single(alpha_y > k)
    disp('Stable in x')
elseif single(alpha_x > k) && single(alpha_y < k)
    disp('Stable in y')
end

% disp(single(k))

toc

clearvars