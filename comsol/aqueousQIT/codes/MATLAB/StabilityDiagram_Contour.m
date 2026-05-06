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
% V_AC = 10 / 2; % [Vpeak]; 2V_AC = AC Voltage [Vpeak]; 4V_AC = AC Voltage [Vpeak-to-peak]
% U_DC = 1e-6 / 2; % [V]
% r_0 = 1e-6; % [m]
% f = 1e9; % [Hz] NOTE: decreasing magnitude by one AND increasing magnitude by one both result in instability; there is an optimal f for a given k
% atomic_numbers = [3];
% [q,a,k,omega,Q_zeta_ratio] = Parameters(V_AC,U_DC,r_0,f,atomic_numbers);
% q_parametertester = q;
% a_parametertester = a;
% clear('q','a')
% q_max = q_parametertester * 1.25;
% q_min = 0;
% a_min = -a_parametertester * 1.25;
% a_max = -a_min;
% q_increment = (q_max - q_min) / 10;
% a_increment = (a_max - a_min) / 10;
% a_constraint = V_AC / U_DC;

% % Na+
% V_AC = 24 / 2; % [Vpeak]; 2V_AC = AC Voltage [Vpeak]; 4V_AC = AC Voltage [Vpeak-to-peak]
% U_DC = 2e-6 / 2; % [V]
% r_0 = 1e-3; % [m]
% f = 1e9; % [Hz] NOTE: decreasing magnitude by one AND increasing magnitude by one both result in instability; there is an optimal f for a given k
% atomic_numbers = [11];
% [q,a,k,omega,Q_zeta_ratio] = Parameters(V_AC,U_DC,r_0,f,atomic_numbers);
% q_parametertester = q;
% a_parametertester = a;
% clear('q','a')
% q_max = q_parametertester * 1.25;
% q_min = 0;
% a_min = -a_parametertester * 1.25;
% a_max = -a_min;
% q_increment = (q_max - q_min) / 100;
% a_increment = (a_max - a_min) / 100;
% a_constraint = V_AC / U_DC;

% % 2micron microsphere
% V_AC = 40 / 2; % [Vpeak]; 2V_AC = AC Voltage [Vpeak]; 4V_AC = AC Voltage [Vpeak-to-peak]
% U_DC = 0.002 / 2; % [V]
% r_0 = 5.531e-3; % [m]
% f = 1e4; % [Hz] NOTE: decreasing magnitude by one AND increasing magnitude by one both result in instability; there is an optimal f for a given k
% atomic_numbers = [124]; % 124-128
% [q,a,k,omega,Q_zeta_ratio] = Parameters(V_AC,U_DC,r_0,f,atomic_numbers);
% q_parametertester = q;
% a_parametertester = a;
% clear('q','a')
% q_max = q_parametertester * 1.25;
% q_min = 0;
% a_min = -a_parametertester * 1.25;
% a_max = -a_min;
% q_increment = (q_max - q_min) / 100;
% a_increment = (a_max - a_min) / 100;
% a_constraint = V_AC / U_DC;

% 10micron microsphere
V_AC = 10 / 2; % [Vpeak]; 2V_AC = AC Voltage [Vpeak]; 4V_AC = AC Voltage [Vpeak-to-peak]
U_DC = 0.002 / 2; % [V]
r_0 = 5.531e-3; % [m]
f = 1e4; % [Hz] NOTE: decreasing magnitude by one AND increasing magnitude by one both result in instability; there is an optimal f for a given k
atomic_numbers = [121]; % 119-123
[q,a,k,omega,Q_zeta_ratio] = Parameters(V_AC,U_DC,r_0,f,atomic_numbers);
q_parametertester = q;
a_parametertester = a;
clear('q','a')
q_max = q_parametertester * 1.25;
q_min = 0;
a_min = -a_parametertester * 1.25;
a_max = -a_min;
q_increment = (q_max - q_min) / 100;
a_increment = (a_max - a_min) / 100;
a_constraint = V_AC / U_DC;

%% Other initializations

% Determine how many elements are in the q-vector based on the step size (increment) and minimum and maximum values
n_q = round((q_max - q_min) / q_increment + 1);

% Establish the q vector
q = q_min:q_increment:q_max;

% Determine how many elements are in the a-vector based on the step size (increment) and minimum and maximum values
n_a = round((a_max - a_min) / a_increment + 1);
% Establish the a-vector
a = a_min:a_increment:a_max;

% a_points = -q ./ a_constraint;

%% Calculate and plot curves

if k == 0
    % Initialize an arbitrary accuracy limit parameter
    relative_error = 1e-7; % [-]

    % Initialize the n-vector
    n = [0];
    % Determine how many elements are in the n-vector
    n_points = length(n);     

    % a2n

    % Initialize the a2n array
    a2n = NaN(n_points,n_q);
    % Calculate the first column of the a2n array
    a2n(:,1) = (2 .* n(:)).^2;

    for n_index = 1:n_points
        clear a2n_approx
        for q_index = 2:n_q
            convergence1 = 0;
            u = 1;    
            a2n_approx(1,1) = a2n(n_index,q_index - 1);
            while convergence1 == 0
                convergence2 = 0;
                s = 1;
                r_a2n(1) = 2;
                while convergence2 == 0   
                    d_a2n(r_a2n(s)+1) = 0;
                    for t = r_a2n(s):-1:2
                        d_a2n(t) = -(1/(0.5 * k^2 + 4 * r_a2n(t-1)^2)) * q(q_index) / (1 - (1/(0.5 * k^2 + 4 * r_a2n(t-1)^2)) * (a2n_approx(u,1) - q(q_index) * d_a2n(t+1)));
                        fprintf('%16.15f\n',d_a2n(t))
                    end 
                    a2n_approx(u,s+1) = -2 / (0.5 * k^2 + 4) * q(q_index)^2 / (1 - 1 / (0.5 * k^2 + 4) * (a2n_approx(u,1) - q(q_index) * d_a2n(2)));
                    fprintf('%16.15f\n',a2n_approx(u,s+1))
                    if abs(a2n_approx(u,s+1) - a2n_approx(u,s)) / abs(a2n_approx(u,s+1)) <= relative_error
                        a2n_approx(u+1,1) = (a2n_approx(u,s+1) + a2n_approx(u,1)) / 2; 
                        fprintf('%16.15f\n',a2n_approx(u+1,1))
                        convergence2 = 1;
                    else
                    r_a2n(s+1) = r_a2n(s) + 1;
                    s = s + 1; 
                    end                        
                end
                if abs(a2n_approx(u+1,1) - a2n_approx(u,1)) / abs(a2n_approx(u+1,1)) <= relative_error
                    a2n(n_index,q_index) = a2n_approx(u+1,1);
                    convergence1 = 1;
                else
                end 
                u = u + 1;
            end
        end
    end

    % b2n+1

    % Initialize the b_2nplus1 array
    b2nplus1 = NaN(n_points,n_q);
    % Calculate the first column of the b_2n+1 array
    b2nplus1(:,1) = (2 .* n(:) + 1).^2;

    for n_index = 1:n_points
        clear b2nplus1_approx
        for q_index = 2:n_q
            convergence1 = 0;
            u = 1;    
            b2nplus1_approx(1,1) = b2nplus1(n_index,q_index - 1);
            while convergence1 == 0
                convergence2 = 0;
                s = 1;
                r_b2nplus1(1) = 1;    
                while convergence2 == 0     
                    d_b2nplus1(r_b2nplus1(s)+2) = 0;
                    for v = r_b2nplus1(s):-1:1
                        d_b2nplus1(v) = q(q_index) / ((b2nplus1_approx(u,1) - (2 * r_b2nplus1(v) + 1)^2 - 0.5 * k^2 - q(q_index) * d_b2nplus1(v+2)));
                        fprintf('%16.15f\n',d_b2nplus1(v))
                    end 
                    b2nplus1_approx(u,s+1) = q(q_index) * d_b2nplus1(1) - q(q_index) + 1 + 0.5 * k^2;
                    fprintf('%16.15f\n',b2nplus1_approx(u,s+1))
                    if abs(b2nplus1_approx(u,s+1) - b2nplus1_approx(u,s)) / abs(b2nplus1_approx(u,s+1)) <= relative_error
                        b2nplus1_approx(u+1,1) = (b2nplus1_approx(u,s+1) + b2nplus1_approx(u,1)) / 2; 
                        fprintf('%16.15f\n',b2nplus1_approx(u+1,1))
                        convergence2 = 1;
                    else
                    r_b2nplus1(s+1) = r_b2nplus1(s) + 1;
                    s = s + 1; 
                    end                             
                end
                if abs(b2nplus1_approx(u+1,1) - b2nplus1_approx(u,1)) / abs(b2nplus1_approx(u+1,1)) <= relative_error
                    b2nplus1(n_index,q_index) = b2nplus1_approx(u+1,1);
                    convergence1 = 1;
                else
                end 
                u = u + 1;
            end
        end
    end

    % Plot        
    hold on
    figure(1)
    plot(q,a2n(1,:),'w',q,-a2n(1,:),'w',q,b2nplus1(1,:),'w',q,-b2nplus1(1,:),'w')        
    % xlabel('q')
    xlim([0 1])
    % ylabel('a','Rotation',0)
    ylim([-0.25 0.25])
    % title('Stability Diagram, k = 0')
    hold on
    q_shade = [q,fliplr(q)];
    x_area = [a2n(1,:),fliplr(b2nplus1(1,:))];
    fill(q_shade,x_area,[1 0 0],'EdgeColor','r','EdgeAlpha',0.6,'FaceAlpha',0.25)
    patch([q_min q_max q_max q_min],[a_min a_min a_max a_max],'w','FaceAlpha',0.4,'EdgeColor','k','EdgeAlpha',0)
    y_area = [-a2n(1,:),fliplr(-b2nplus1(1,:))];
    fill(q_shade,y_area,[0 0 1],'EdgeColor','b','FaceAlpha',0.4)
    hold off

else

    r_max = 100;
    r_min = -r_max;
    r_points = r_max - r_min + 1;
    r = linspace(r_min,r_max,r_points);

    % Initialize the infinite determinant matrix
    inf_matrix_zero = eye(r_points);

    for n = 1:n_a
        for m = 1:n_q

            % Calculate all of the xi's necessary for the infinite determinant matrix
            xi = q(m) ./ (4 * r.^2 - a(n) + k^2 / 4);

            for f = 1:r_points - 1
                y(f) = f + 1;
                x(f) = f - 1;

                % Left diagonal of infinite matrix                
                inf_matrix_zero(y(f),f) = xi(y(f));

                % Right diagonal of infinite matrix
                inf_matrix_zero(f,y(f)) = xi(f);
            end

            % Calculate the determinant of the infinite matrix and its derivative
            inf_det_zero = det(inf_matrix_zero);

            % Calculate alpha

            if single(k) > 380
                alpha(n,m) = 4 / pi * asinh(sqrt(-inf_det_zero * sin(sym(pi / 2 * sqrt(a(n) - k^2 / 4)))^2));
            else
                alpha(n,m) = 4 / pi * asinh(sqrt(-inf_det_zero .* sin(pi / 2 * sqrt(a(n) - k^2 / 4))^2));
            end  

            % disp(m)
            % toc

        end

        % disp(n)
        % toc

    end

    alpha = real(alpha);
    
    hold on    

    %% Plot
    figure(1)
    xlabel('q')
    ylabel('a','Rotation',0)
    xlim([q_min q_max])
    ylim([a_min a_max])
    title({"Stability Diagram, k = " + k + ""})
    % title({"Stability Diagram, k = " + k + " (Maximum Size, Maximum Charge)"})
    % title("Stability Diagram Comparison (2 vs. 0.5 micron microspheres)")
    % title("Stability Diagram Comparison (Li^+ and Na^+)")
    ax1 = axes;
    [Cx,hx] = contourf(q,a,alpha,[-inf k inf],'ShowText',false,'EdgeColor','r','FaceAlpha',0.4); % [-inf k inf] is used for shading purposes. [k k] provides lines only
    patch([q_min q_max q_max q_min],[a_min a_min a_max a_max],'w','FaceAlpha',0.4)
    view(2)
    ax2 = axes;
    [Cy,hy] = contourf(q,-a,alpha,[-inf k inf],'ShowText',false,'EdgeColor','b','FaceAlpha',0.4);
    linkaxes([ax1,ax2])
    ax1.Visible = 'off';
    ax2.Visible = 'off';
    cmapx = [1 0 0 ; 1 1 1];
    cmapy = [0 0 1 ; 1 1 1];
    colormap(ax1,cmapx)
    colormap(ax2,cmapy)
    patch([q_min q_max q_max q_min],[a_min a_min a_max a_max],'k','FaceAlpha',0) 

    ax1 = axes;
    [Cx,hx] = contourf(q,a,alpha,[-inf k inf],'ShowText',false,'EdgeColor','r','FaceAlpha',0);
    patch([q_min q_max q_max q_min],[a_min a_min a_max a_max],'w','FaceAlpha',0)
    view(2)
    ax2 = axes;
    [Cy,hy] = contourf(q,-a,alpha,[-inf k inf],'ShowText',false,'EdgeColor','b','FaceAlpha',0);
    linkaxes([ax1,ax2])
    ax1.Visible = 'off';
    ax2.Visible = 'off';
    cmapx = [1 0 0 ; 1 1 1];
    cmapy = [0 0 1 ; 1 1 1];
    colormap(ax1,cmapx)
    colormap(ax2,cmapy)
    hold off
    hold on
    text(q_parametertester,a_parametertester * 0.95,"(" + sprintf('%.2e',q_parametertester) + "," + sprintf('%.2e',a_parametertester) + ")",'HorizontalAlignment','center','VerticalAlignment','top','Interpreter','latex','FontSize',16,'Color','k')
    point = plot(q_parametertester,a_parametertester,'.k','MarkerSize',15);

    linkaxes([ax1,ax2])
    hold off

end

toc