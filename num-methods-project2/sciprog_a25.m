% kundyz muktar, matlab r2023a, sci prog a2, ex 5
clear; clc; close all;

function [z, zvec] = newton(f_handle, df_handle, z0)
    % implements newton's method for f_handle: Rd to Rd w jacobian df_handle
    % columns of zvec are the iterates. z is the final one.
    
    abs_tol = 1e-12;
    max_iters = 1000;
    d = length(z0);
    zvec = zeros(d, max_iters + 1);
    zvec(:, 1) = z0;

    for iter_index = 1:max_iters
        current_z = zvec(:, iter_index);
        % newton step
        zvec(:, iter_index + 1) = current_z - df_handle(current_z)\f_handle(current_z);
        if norm(zvec(:, iter_index + 1) - current_z) <= abs_tol
            break;
        end
    end

    zvec = zvec(:, 1:iter_index + 1); % when method stops early, reset the zvec size
    z = zvec(:, end); % extracting the final sol
end
