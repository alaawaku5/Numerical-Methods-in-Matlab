clc;    
clear;  


f1 = @(t) 0;        % Constant forcing function
f2 = @(t) -1.0;       % Constant coefficient for y
y0 = 1; % iv
exact_sol= @(t)exp(-t);


%______________________________
L =10; 

error_explicit = zeros(1, L);
error_implicit = zeros(1, L);
error_midpoint = zeros(1, L);
h_vals = zeros(1, L);

for i=1:L
    
    h = 2^-i;
    h_vals(i)=h; % store vals of h for eaech L for future use
    tfirst=0; tlast=1; 
    t = tfirst:h:tlast; % % time steps from 0 to 1 in increments of h
   
    y_approx_explicit = explicit_Euler(t, f1, f2, y0);
    y_approx_implicit = implicit_Euler(t, f1, f2, y0);
    y_approx_midpoint = implicit_Midpoint(t, f1, f2, y0);
    y_exact = exact_sol(t);
    
    error_explicit(i) = max(abs(y_approx_explicit - y_exact));
    error_implicit(i) = max(abs(y_approx_implicit - y_exact));
    error_midpoint(i) = max(abs(y_approx_midpoint - y_exact));
end;

errorMatrix = [error_explicit; error_implicit; error_midpoint];
disp(h_vals)
disp(error_explicit)
%______________________________

figure;
loglog(h_vals, error_explicit, 'm', 'DisplayName','Explicit Euler');
hold on;
loglog(h_vals, error_implicit, 'g', 'DisplayName','Implicit Euler');
loglog(h_vals, error_midpoint, 'r', 'DisplayName','Implicit Midpoint');
hold off;
xlabel('Step size h');
ylabel('Max Error');
title('Part (b): Error vs. Step Size (Log-Log Plot)');
legend('Location','best');
grid on;
%______________________________
orderMatrix = zeros(3, L-1);
for i = 1:3
    for j = 1:(L-1)
        orderMatrix(i, j) = log(errorMatrix(i,j) / errorMatrix(i,j+1)) / log(2);
    end
end
%______________________________
fprintf('Error matrix (rows: Exp, Imp, Mid):\n');
disp(errorMatrix);
fprintf('Order of convergence matrix (rows: Exp, Imp, Mid):\n');
disp(orderMatrix);

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
