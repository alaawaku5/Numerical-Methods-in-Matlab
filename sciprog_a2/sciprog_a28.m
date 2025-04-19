% kundyz muktar, matlab r2023a, sci prog a2, ex 8

clear; clc; close all;
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
