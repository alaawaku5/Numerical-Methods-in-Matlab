% kundyz muktar, matlab r2023a, sci prog a2, ex 3

clear; clc;

f = @(t, y) [-5*y(1);  % example: y1' = -2*y1
     y(2)^2 - y(1)+16  % example: y2' = y1^2 - y2
     ];

T = 1.0;     
N = 10;      

y0 = [1; 0];

[timeExp, y_explicit] = explicit_euler(f, y0, T, N);
[timeImp, y_implicit] = implicit_euler(f, y0, T, N);
[timeMid, y_midpoint] = implicit_midpoint(f, y0, T, N);

% Display results (you can also plot them)
disp('Explicit Euler (y_explicit):');
disp(y_explicit);
disp('Implicit Euler (y_implicit):');
disp(y_implicit)
disp('Implicit Midpoint (y_midpoint):');
disp(y_midpoint);
%% 
% basically nothing changes
function [t_values, y_values] = explicit_euler(f, y0, T, N)
    h = (T - 0)/N;
    d = length(y0);
    y_values = zeros(d, N+1);
    t_values = linspace(0, T, N+1);
    y_values(:, 1) = y0;
    
    for j = 1:N
        t_j = t_values(j);
        y_values(:, j+1) = y_values(:, j) + h * f(t_j, y_values(:, j));
    end
end

%% 
% we change the line to find z, we call fixed point
% previously, we could find y_k+1 exactly for nonlinear odes', which gave
% us immediate convergence
function [t_values, y_values] = implicit_euler(f, y0, T, N)
    h = (T - 0)/N;
    d = length(y0);
    y_values = zeros(d, N+1);
    t_values = linspace(0, T, N+1);
    y_values(:, 1) = y0;
    
    for j = 1:N
        t_j = t_values(j);
        y_j = y_values(:, j);
        t_next = t_j + h; % next time step
        %%%%%%%%%%%%%%%%%%
        G = @(z) y_j + h * f(t_next, z);  % define G as given for fixed point 
        [z, ~] = fixedpoint(G, y_j); % solve z = G(z) with z0 = y_j
        %%%%%%%%%%%%%%%%%%
        y_values(:, j+1) = z; %last approx
    end
end

%% 
% we change the line to find z, we call fixed point
function [t_values, y_values] = implicit_midpoint(f, y0, T, N)
    h = (T - 0)/N;
    d = length(y0);
    y_values = zeros(d, N+1);
    t_values = linspace(0, T, N+1);
    y_values(:, 1) = y0;
    
    for j = 1:N
        t_j = t_values(j);
        y_j = y_values(:, j);
        
        % define G for midpoint: z -> y_j + h*f(t_j + h/2, (y_j + z)/2)
        G = @(z) y_j + h * f(t_j + h/2, (y_j + z)/2);
        
        % solve via fixedpoint with initial guess y_j
        [z, ~] = fixedpoint(G, y_j);
        
        y_values(:, j+1) = z;
    end
end
%% 

% fixed point (from Exercise 1)
function [z, zvec] = fixedpoint(G, z0)

    tau_abs = 1e-12;
    Jmax = 1000; 

    % dimension of z0
    d = length(z0);

    zvec = zeros(d, Jmax + 1);
    zvec(:, 1) = z0;

    for j = 1:Jmax
        zvec(:, j + 1) = G(zvec(:, j));
        if norm(zvec(:, j + 1) - zvec(:, j)) <= tau_abs
            break;
        end
    end

    zvec = zvec(:, 1:j + 1); % trim unused columns
    z = zvec(:, end); % final approx
end


