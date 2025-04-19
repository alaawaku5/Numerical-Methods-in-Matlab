% kundyz muktar, matlab r2023a, sci prog a2, ex 2
clear; clc;

% from problem def
a = 0.5;
G = @(x) exp(-a * abs(x));
x0 = 100;


[z, zvec] = fixedpoint(G, x0);

disp('Approximate solution x* such that G(x*) = x*:');
disp(z);

% show how many iterations we used
disp('Number of iterations:');
disp(size(zvec, 2));


% disp('All iterates (zvec):');
% disp(zvec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%same as ex1
function [z, zvec] = fixedpoint(G, z0)
    tau_abs = 1e-12;
    Jmax = 1000; 

    zvec = zeros(1, Jmax + 1); % to store all iteates
    zvec(:, 1) = z0; % first col is the initial guess

    for j = 1:Jmax
        zvec(:, j + 1) = G(zvec(:, j));
        if norm(zvec(:, j + 1) - zvec(:, j)) <= tau_abs
            break;
        end
    end

    zvec = zvec(:, 1:j + 1); % trim unused columns
    z = zvec(:, end); % final approx
end
