clc;  

N = 10;    % Maximum exponent for h = 2^(-n)
x = 0;     % Evaluation point

% Define the functions and their names in cell arrays
funcs = { @(x) exp(x), @(x) sin(x), @(x) x.^(3/2) };
fnames = { 'exp(x)', 'sin(x)', 'x^(3/2)' };
fprime_zero = [1,1,0];

errors_matrix = zeros(3, N+1);
order_matrix = zeros(3,N);
figure; hold on;  % PLOT
% computing one-sided difference quotient
for i = 1:3
    
    fprintf('Function: %s\n', fnames{i});

    [h_values, dq_values] = oneSidedDQ(funcs{i}, x, N);

    % ERROR 
    error = abs(dq_values - fprime_zero(i));
    errors_matrix(i,:) = error;
    fprintf('\n Error vector for %s:\n', fnames{i}); %display error
    error
    
    % ORDER P
    for j = 1:N  % Loop through columns (excluding the last)
        order_matrix(i,j) = (log(error(j))-log(error(j+1)))/log(2);
    end


    loglog(h_values, error, '-o', 'LineWidth', 1.5, 'DisplayName', fnames{i}); % PLOT
    fprintf("------------------------------------------\n");
end

% Formatting the plot
xlabel('Step size h');
ylabel('Error e_h');
title('Error vs. Step Size (Log-Log Plot)');
legend show;
grid on;
hold off;


fprintf('\nError Table:\n');
disp(errors_matrix); 
fprintf("------------------------------------------\n");
fprintf('\nOrder of convergence matrix:\n');
disp(order_matrix);






%% Local Function Definition
function [h_vals, dq_vals] = oneSidedDQ(f, x, N)

% allocating space
    h_vals = zeros(1, N+1);
    dq_vals = zeros(1, N+1);

     % compute diff quotient for each n in N
    for n = 0:N
        h = 2^(-n);
        h_vals(n+1) = h; % storing vals in the array for future retrieval
        dq_vals(n+1) = (f(x + h)-f(x))/ h; % arrays in matlab start from 1!!!
    end 

    fprintf('   n       h_values                  dq_values\n');
    for n = 0:N
        fprintf('%3d   %12.5g      %14.8g\n', n, h_vals(n+1), dq_vals(n+1));
    end
end
