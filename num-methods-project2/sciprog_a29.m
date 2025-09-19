% kundyz muktar, matlab r2023a, sci prog a2, ex 9
clear; clc; close all;

% Define the ODE so that the exact solution is [cos(t); sin(t)]
f = @(t, y) [y(2); -y(1)];
exact_sol = @(t) [cos(t); sin(t)];

t0 = 0;
T = 2*pi;
y0 = [1; 0];

n_values = 1:1:25;
% disp(n_values)

errors_uniform = zeros(size(n_values));

% Uniform Euler method
for i = 1:length(n_values)
    [t_vals_uni, y_vals_uni] = uniform_euler(f, t0, y0, T, n_values(i));

    final_approx_uni = y_vals_uni(:, end);
    exact_at_end = exact_sol(T);

    errors_uniform(i) = norm(final_approx_uni - exact_at_end, inf);
end

% Adaptive Euler with local error tolerance ~ 1/n_values, starting with h = 1
errors_adaptive = zeros(size(n_values));
steps_adaptive = zeros(size(n_values));
h_init = 5.0;

for i = 1:length(n_values)
    N = n_values(i);
    tolerance = 2^(-N);
    [t_vals_adapt, y_vals_adapt] = adaptiveTimeStepping(f, t0, y0, T, h_init, tolerance, 'explicit_euler');

    final_approx_adapt = y_vals_adapt(:, end);
    exact_at_end = exact_sol(T);

    errors_adaptive(i) = norm(final_approx_adapt - exact_at_end, inf);
    steps_adaptive(i) = length(t_vals_adapt) - 1; % number of steps used
end

% Plot errors vs. number of steps in a loglog scale
figure;
loglog(n_values, errors_uniform, 'o-'); hold on;
loglog(steps_adaptive, errors_adaptive, 's-');
xlabel('Number of Steps');
ylabel('Error at T = 2\pi');
legend('Uniform Euler', 'Adaptive Euler', 'location', 'best');
title('Error vs. Number of Steps (loglog)');
grid on;


%% 

function [t0, y] = uniform_euler(f, t0, y0, T, N)
    % does explicit euler w/ uniform step
    % h = 2^(-N);
    h = (T-t0)/N;
    t0 = linspace(t0, T, N + 1);
    d = length(y0);
    y = zeros(d, N + 1);
    y(:, 1) = y0;
    for j = 1:N
        y(:, j+1) = y(:, j) + h * f(t0(j), y(:, j));
    end
end

%%
% kundyz muktar, matlab r2023a, sci prog a2, ex 8

function [t, y] = adaptiveTimeStepping(f, t0, y0, T, h0, tau, solveODE)
    t_curr = t0;
    y_curr = y0;
    h_curr = h0;

    t = t_curr;
    y = y_curr;

    while t_curr < T
        if t_curr + h_curr > T
            h_curr = T - t_curr; %adjust h such that it reaches T exactly
        end
        
        y_big = single_step(f, t_curr, y_curr, 2*h_curr, solveODE); % one big step
        
        y_half  = single_step(f, t_curr,        y_curr,  h_curr, solveODE); % two half steps
        y_half2 = single_step(f, t_curr + h_curr, y_half, h_curr, solveODE);

        err_est = norm(y_big - y_half2, inf);  

        if err_est <= tau % if the desired accuracy was achieved
            t_curr = t_curr + 2*h_curr;
            y_curr = y_half2;
            
            t(end+1) = t_curr;
            y(:,end+1) = y_curr;
        else
            h_curr = 0.5 * h_curr; %halve the step size and repeat
        end
    end
end





%%
%choosing the solver method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y_next = single_step(f, t_now, y_now, h, solveODE)
    % does one step of chosen solveODE
    switch lower(solveODE)
        case 'explicit_euler'
            y_next = y_now + h * f(t_now, y_now);
        case 'implicit_midpoint'
            y_next = implicit_midpoint_solve(f, t_now, y_now, h);
        otherwise
            error('unknown solveODE');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y_sol = implicit_midpoint_solve(f, t_now, y_now, h)
    % solves y = y_now + h*f(t_now + h/2, (y_now + y)/2) by fixedpoint
    abs_tol = 1e-12;
    max_iters = 1000;
    d = length(y_now);
    guess_vals = zeros(d, max_iters+1);
    guess_vals(:,1) = y_now;
    for i = 1:max_iters
        guess_vals(:,i+1) = y_now + h * f(t_now + h/2, (y_now + guess_vals(:,i))/2);
        if norm(guess_vals(:,i+1) - guess_vals(:,i)) <= abs_tol
            break;
        end
    end
    y_sol = guess_vals(:, i+1);
end
