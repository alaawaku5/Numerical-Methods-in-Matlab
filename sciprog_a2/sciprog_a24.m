% kundyz muktar, matlab r2023a, sci prog a2, ex 4
clear; clc; close all;

% set params
k = 1;
y0 = 0.5;
T = 5;

f = @(t,y) k * y .* (1 - y);
exact = @(t) 1./(1 + ((1/y0) - 1) .* exp(-k .* t));
N      = 1:4;
h_vals = 2.^(-N);

%%
% storing errors per h val
err_exp = zeros(size(h_vals));
err_imp = zeros(size(h_vals));
err_mid = zeros(size(h_vals));

%%
% explicit euler error vs h
figure; hold on
for i = 1:length(h_vals)
    h = h_vals(i);
    N = round(T/h);
    [t,y] = explicit_euler(f, y0, T, N);
    err_exp(i) = max(abs(y - exact(t)));
    disp(err_exp(i));
    plot(h, err_exp(i), '-.', 'displayname', sprintf('h=%.3f', h));
end
xlabel('h'), ylabel('error'), title('explicit euler'), legend, grid on

%%
% implicit euler iters vs time
figure; hold on
for i = 1:length(h_vals)
    h = h_vals(i);
    N = round(T/h);
    [t,y,iters] = implicit_euler(f, y0, T, N);
    err_imp(i) = max(abs(y - exact(t)));
    plot(t(2:end), iters, '-.', 'displayname', sprintf('h=%.3f', h));
end
xlabel('t'), ylabel('iters'), title('implicit euler'), legend, grid on

%%
% implicit midpoint iters vs time
figure; hold on
for i = 1:length(h_vals)
    h = h_vals(i);
    N = round(T/h);
    [t,y,iters] = implicit_midpoint(f, y0, T, N);
    err_mid(i) = max(abs(y - exact(t)));
    plot(t(2:end), iters, '-.', 'displayname', sprintf('h=%.3f', h));
end
xlabel('t'), ylabel('iters'), title('implicit midpoint'), legend, grid on

%%
% compare errors vs step size
figure;
loglog(h_vals, err_exp, 's-', h_vals, err_imp, 's-', h_vals, err_mid, 's-');
legend('Explicit Euler','Implicit Euler','Implicit Midpoint','Location','best');
xlabel('Step Size');
ylabel('Max Error');
title('Error Comparison for Logistic ODE');
grid on;









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% from ex 3
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
function [t_values, y_values, iter_count] = implicit_euler(f, y0, T, N)
    h = (T - 0)/N;
    d = length(y0);
    y_values = zeros(d, N+1);
    t_values = linspace(0, T, N+1);
    iter_count = zeros(1, N);
    y_values(:, 1) = y0;
    
    for j = 1:N
        t_j = t_values(j);
        y_j = y_values(:, j);
        t_next = t_j + h; % next time step
        %%%%%%%%%%%%%%%%%%
        G = @(z) y_j + h * f(t_next, z);  % define G as given for fixed point 
        [z, iters] = fixedpoint(G, y_j); % solve z = G(z) with z0 = y_j
        % if j==1
        %     [dim, len] = size(iters);
        %     disp(len)
        % 
        % end
        [~, len] = size(iters);
        iter_count(j) = len;
        %%%%%%%%%%%%%%%%%%
        y_values(:, j+1) = z; %last approx
    end
end

%%
% we change the line to find z, we call fixed point
function [t_values, y_values, iter_count] = implicit_midpoint(f, y0, T, N)
    h = (T - 0)/N;
    d = length(y0);
    y_values = zeros(d, N+1);
    t_values = linspace(0, T, N+1);
    iter_count = zeros(1, N);
    y_values(:, 1) = y0;
    
    for j = 1:N
        t_j = t_values(j);
        y_j = y_values(:, j);
        
        % define G for midpoint: z -> y_j + h*f(t_j + h/2, (y_j + z)/2)
        G = @(z) y_j + h * f(t_j + h/2, (y_j + z)/2);
        
        % solve via fixedpoint with initial guess y_j
        [z, iters] = fixedpoint(G, y_j);

        [~, len] = size(iters);
        iter_count(j) = len;
        
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


