function [q,a,k,omega,Q_zeta_ratio] = Parameters(V_AC,U_DC,r_0,f,atomic_numbers)

% Import atomic and ionic data
data = readtable("ElementandMaterialDatabase.xlsx");

% Extract the element names
ion_name = data{:,"Element"};

% Extract the ionic charges
ionic_Q = data{:,"Charge"};

% Extract the ionic masses
ionic_M = data{:,"IonicMass"};

% Extract the ionic radii
ionic_r = data{:,"IonicRadius"};

% Extract the zeta values
zeta = data{:,"Zeta"};

% Angular frequency
omega = 2 * pi * f; % [Hz]

% Dynamic viscosity of water
eta = 0.001; % [N*s/m^2]

for n = 1:length(atomic_numbers)
    a(n) = 8 * ionic_Q(atomic_numbers(n)) * U_DC / (ionic_M(atomic_numbers(n)) * r_0^2 * omega^2);
    q(n) = 4 * ionic_Q(atomic_numbers(n)) * V_AC / (ionic_M(atomic_numbers(n)) * r_0^2 * omega^2);
    zeta(n) = zeta(atomic_numbers(n));
    Q_zeta_ratio(n) = ionic_Q(atomic_numbers(n)) / zeta(n);
    k(n) = 2 * zeta(n) / (ionic_M(atomic_numbers(n)) * omega);
    if k(n) > 380
        clear k
        k(n) = sym(2 * zeta(n) / (ionic_M(atomic_numbers(n)) * omega));
    end
    name(n) = ion_name(atomic_numbers(n)); 
end

format('long')

T = table(transpose(a),transpose(q),transpose(k),'VariableNames',["a","q","k"],'RowNames',name);
disp(T)

end

