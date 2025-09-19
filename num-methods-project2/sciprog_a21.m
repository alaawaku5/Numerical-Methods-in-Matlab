% kundyz muktar, matlab r2023a, sci prog a2, ex 1
clear; clc;

G = @(x) [cos(x(1)); x(2) / (2 + x(1)^2)]; % ex
z0 = [1; 1]; % initial guess

[z, zvec] = fixedpoint(G, z0);


disp('final z (approx fixed point):');
disp(z);

disp('all iterates:');
disp(zvec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z, zvec] = fixedpoint(G, z0)
    tau_abs = 1e-12;
    Jmax = 1000;

    d = length(z0);
    zvec = zeros(d, Jmax + 1);
    zvec(:, 1) = z0;

    for j = 1:Jmax
        zvec(:, j + 1) = G(zvec(:, j));
        if norm(zvec(:, j + 1) - zvec(:, j)) <= tau_abs
            break; %break if the change is too small
        end
    end

    zvec = zvec(:, 1:j + 1); % trim unused columns
    z = zvec(:, end); % final approx
end
