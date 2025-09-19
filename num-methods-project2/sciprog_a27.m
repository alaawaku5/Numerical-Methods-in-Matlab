% kundyz muktar, matlab r2023a, sci prog a2, ex 7
% apply one step methods with a. fixed point b. newton's to lotka volterra
% similiar as problem 3, but instead of fixed point we use newton

clear; clc; close all;

T = 15;
N = 150;
step_size = T / N;


f = @(y) [y(1)*(1-y(2)); -y(2)*(1-y(1))];

f_jac = @(y) [ 1-y(2),  -y(1);    y(2),    y(1)-1 ];


y0 = [1; 0.5];
%%
% (a) exercise 3 methods w/ fixedpoint
[t_exp, y_e] = explicit_euler(f, y0, T, N);
[t_ie, y_ie, iter_count_ie] = implicit_euler(f, y0, T, N);
[t_im, y_im, iter_count_im] = implicit_midpoint(f, y0, T, N);

% (b) newton-based methods
[t_ien, y_ien, iter_count_ien] = implicit_euler_newton(f, f_jac, y0, T, N);
[t_imn, y_imn, iter_count_imn] = implicit_midpoint_newton(f, f_jac, y0, T, N);

%%
% compare iteration counts for the implicit methods
figure;
plot(t_ie(2:end), iter_count_ie); hold on;
plot(t_ien(2:end), iter_count_ien);
xlabel('time'); ylabel('iterations');
legend('implicit euler (fixedpoint)','implicit euler (newton)','location','best');
title('compare iteration counts: euler');

figure;
plot(t_im(2:end), iter_count_im); hold on;
plot(t_imn(2:end), iter_count_imn);
xlabel('time'); ylabel('iterations');
legend('implicit midpoint (fixedpoint)','implicit midpoint (newton)','location','best');
title('compare iteration counts: midpoint');



% from previous exercises
%% explicit euler
function [t_values, y_values] = explicit_euler(f_handle, y0, T, N)
    h = T / N;
    t_values = linspace(0, T, N + 1);
    d = length(y0);
    y_values = zeros(d, N + 1);
    y_values(:,1) = y0;
    for j = 1:N
        y_values(:, j+1) = y_values(:, j) + h * f_handle(y_values(:, j));
    end
end

%% implicit euler (fixedpoint)
function [t_values, y_values, iter_count] = implicit_euler(f_handle, y0, T, N)
    h = T / N;
    t_values = linspace(0, T, N + 1);
    d = length(y0);
    y_values = zeros(d, N + 1);
    y_values(:,1) = y0;
    iter_count = zeros(1, N);
    for j = 1:N
        y_j = y_values(:, j);
        g_func = @(z) y_j + h * f_handle(z);
        [z_val, ~, iter_used] = fixedpoint(g_func, y_j);
        y_values(:, j+1) = z_val;
        iter_count(j) = iter_used;
    end
end

%% implicit midpoint (fixedpoint)
function [t_values, y_values, iter_count] = implicit_midpoint(f_handle, y0, T, N)
    h = T / N;
    t_values = linspace(0, T, N + 1);
    d = length(y0);
    y_values = zeros(d, N + 1);
    y_values(:,1) = y0;
    iter_count = zeros(1, N);
    for j = 1:N
        y_j = y_values(:, j);
        g_func = @(z) y_j + h * f_handle((y_j + z)/2);
        [z_val, ~, iter_used] = fixedpoint(g_func, y_j);
        y_values(:, j+1) = z_val;
        iter_count(j) = iter_used;
    end
end

%% implicit euler (newton)
function [t_values, y_values, iter_count] = implicit_euler_newton(f_handle, jac_handle, y0, T, N)
    h = T / N;
    t_values = linspace(0, T, N + 1);
    d = length(y0);
    y_values = zeros(d, N + 1);
    y_values(:,1) = y0;
    iter_count = zeros(1, N);
    for j = 1:N
        y_j = y_values(:, j);
        % solve f_e(z)=0 where f_e(z)= z - y_j - h*f_handle(z)
        newton_func = @(z) z - y_j - h * f_handle(z);
        newton_jac = @(z) eye(d) - h * jac_handle(z);
        [z_val, ~, iter_used] = newton(newton_func, newton_jac, y_j);
        y_values(:, j+1) = z_val;
        iter_count(j) = iter_used;
    end
end

%% implicit midpoint (newton)
function [t_values, y_values, iter_count] = implicit_midpoint_newton(f_handle, jac_handle, y0, T, N)
    h = T / N;
    t_values = linspace(0, T, N + 1);
    d = length(y0);
    y_values = zeros(d, N + 1);
    y_values(:,1) = y0;
    iter_count = zeros(1, N);
    for j = 1:N
        y_j = y_values(:, j);
        % solve f_m(z)=0 where f_m(z)= z - y_j - h*f_handle((y_j+z)/2)
        newton_func = @(z) z - y_j - h * f_handle((y_j + z)/2);
        newton_jac = @(z) eye(d) - (h/2) * jac_handle((y_j + z)/2);
        [z_val, ~, iter_used] = newton(newton_func, newton_jac, y_j);
        y_values(:, j+1) = z_val;
        iter_count(j) = iter_used;
    end
end

%% fixedpoint (from ex 1) fixedpoint
function [z_val, z_values, iteration_used] = fixedpoint(g_handle, z0)
    abs_tol = 1e-12;
    max_iters = 1000;
    d = length(z0);
    z_values = zeros(d, max_iters + 1);
    z_values(:,1) = z0;
    for n = 1:max_iters
        z_values(:,n+1) = g_handle(z_values(:,n));
        if norm(z_values(:,n+1) - z_values(:,n)) <= abs_tol
            break;
        end
    end
    iteration_used = n;
    z_values = z_values(:,1:n+1);
    z_val = z_values(:,end);
end

%% newton (from ex 5) Newton's 
function [z_val, z_values, iteration_used] = newton(f_handle, df_handle, z0)
    abs_tol = 1e-12;
    max_iters = 1000;
    d = length(z0);
    z_values = zeros(d, max_iters + 1);
    z_values(:,1) = z0;
    for i = 1:max_iters
        current_z = z_values(:,i);
        z_values(:,i+1) = current_z - df_handle(current_z)\f_handle(current_z);
        if norm(z_values(:,i+1) - current_z) <= abs_tol
            break;
        end
    end
    iteration_used = i;
    z_values = z_values(:,1:i+1);
    z_val = z_values(:,end);
end
