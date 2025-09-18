
clc;    
clear;  

t = 0:0.1:1; % time steps from 0 to 1 in increments of 0.1
f1 = @(t) 2.0;        % Constant forcing function
f2 = @(t) -1.0;       % Constant coefficient for y
y0 = 5; % ic

% Compute the approximate solutions using different methods
y_approx_explicit = explicit_Euler(t, f1, f2, y0);
y_approx_implicit = implicit_Euler(t, f1, f2, y0);
y_approx_midpoint = implicit_Midpoint(t, f1, f2, y0);

% Display the results
disp('Time steps:');
disp(t);

disp('Approximate solution using Explicit Euler:');
disp(y_approx_explicit);

disp('Approximate solution using Implicit Euler:');
disp(y_approx_implicit);

disp('Approximate solution using Implicit Midpoint:');
disp(y_approx_midpoint);

%% Local Function Definitions

function y = explicit_Euler(t, f1, f2, y0)

    L = length(t) - 1;    % number of steps
    y = zeros(1, L+1);    % pre-allocate solution vector
    y(1) = y0;            % initial condition
    
    for n = 1:L
        h = t(n+1) - t(n);        % step size
        fVal = f1(t(n)) + f2(t(n)) * y(n);
        y(n+1) = y(n) + h * fVal;   % Euler update
    end
end

function y = implicit_Euler(t, f1, f2, y0)

    L = length(t) - 1;
    y = zeros(1, L+1);
    y(1) = y0;
    
    for n = 1:L
        h = t(n+1) - t(n);
        % Evaluate f1 and f2 at t(n+1)
        A = f1(t(n+1));
        B = f2(t(n+1));
        y(n+1) = (y(n) + h*A) / (1 - h*B);
    end
end

function y = implicit_Midpoint(t, f1, f2, y0)
    L = length(t) - 1;
    y = zeros(1, L+1);
    y(1) = y0;
    
    for n = 1:L
        h = t(n+1) - t(n);
        tm = (t(n) + t(n+1)) / 2;  % midpoint in time
        % Evaluate f1 and f2 at the midpoint
        A = f1(tm);
        B = f2(tm);
        y(n+1) = (y(n) + h*A + (h*B/2)*y(n)) / (1 - (h*B/2));
    end
end
