clc;
clear;

L = 10;

error_explicit = zeros(1, L);
error_implicit = zeros(1, L);
error_midpoint = zeros(1, L);
h_vals = zeros(1, L);

d = 2;
f1 = @(t) zeros(d,1);
f2 = @(t) [0 1; -1 0];
y0 = [1; 0];

% Exact solution: y(t) = [cos(t); sin(t)]
exact_sol = @(t) [cos(t); sin(t)];

for i = 1:L
    h = 2^(-i); % step size
    h_vals(i) = h; % storing h for future retrieval

    tfirst = 0; tlast = 1; %grid
    t = tfirst:h:tlast;

    % numerical solutiosn
    y_explicit = explicit_Euler(t, f1, f2, y0);
    y_implicit = implicit_Euler(t, f1, f2, y0);
    y_midpoint = implicit_Midpoint(t, f1, f2, y0);

    % Compute exact solution on the same grid
    y_exact = zeros(d, length(t));
    for n = 1:length(t)
        y_exact(:,n) = exact_sol(t(n));
    end

    % error_explicit(i) = max over n of ||y_explicit(:,n) - y_true(:,n)||_
    error_explicit(i) = max(vecnorm(y_explicit - y_exact, inf, 1));
    error_implicit(i) = max(vecnorm(y_implicit - y_exact, inf, 1));
    error_midpoint(i) = max(vecnorm(y_midpoint - y_exact, inf, 1));
end

% Build error matrix and display
errorMatrix = [error_explicit; error_implicit; error_midpoint];
disp('h values:');
disp(h_vals);
disp('Error matrix (rows: exp, imp, mid):');
disp(errorMatrix);

% Plot error vs h in log-log
figure;
loglog(h_vals, error_explicit, 'm', 'DisplayName','Explicit Euler');
hold on;
loglog(h_vals, error_implicit, 'g', 'DisplayName','Implicit Euler');
loglog(h_vals, error_midpoint, 'r', 'DisplayName','Implicit Midpoint');
hold off;
xlabel('Step size h');
ylabel('Max Error (Euclidean norm)');
title('2D System: Error vs. Step Size (Log-Log Plot)');
legend('Location','best');
grid on;

% Compute and display order of convergence
orderMatrix = zeros(3, L-1);
for row = 1:3
    for j = 1:(L-1)
        orderMatrix(row, j) = log(errorMatrix(row,j) / errorMatrix(row,j+1)) / log(2);
    end
end
disp('Order of convergence matrix (rows: Exp, Imp, Mid):');
disp(orderMatrix);

%% Local Function Definitions

function y = explicit_Euler(t, f1, f2, y0)
    d = length(y0);
    L = length(t) - 1;
    y = zeros(d, L+1);
    y(:,1) = y0;

    for n = 1:L
        h = t(n+1) - t(n);
        A = f1(t(n));     % a d-vector (here 0)
        B = f2(t(n));     % a dxd matrix
        y(:,n+1) = y(:,n) + h * (A + B*y(:,n));
    end
end

function y = implicit_Euler(t, f1, f2, y0)
    d = length(y0);
    L = length(t) - 1;
    y = zeros(d, L+1);
    y(:,1) = y0;

    for n = 1:L
        h = t(n+1) - t(n);
        A = f1(t(n+1));
        B = f2(t(n+1));
        I = eye(d);
        % Solve [I - hB]*y(:,n+1) = y(:,n) + h*A
        y(:,n+1) = (I - h*B)\( y(:,n) + h*A );
    end
end

function y = implicit_Midpoint(t, f1, f2, y0)
    d = length(y0);
    L = length(t) - 1;
    y = zeros(d, L+1);
    y(:,1) = y0;

    for n = 1:L
        h = t(n+1) - t(n);
        tm = (t(n) + t(n+1))/2;
        A = f1(tm);
        B = f2(tm);
        I = eye(d);
        % Solve [I - (h/2)*B]*y(:,n+1) = y(:,n) + h*A + (h/2)*B*y(:,n)
        y(:,n+1) = (I - (h/2)*B)\( y(:,n) + h*A + (h/2)*B*y(:,n) );
    end
end
