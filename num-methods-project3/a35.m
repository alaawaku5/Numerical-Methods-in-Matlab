% kundyz muktar, matlab r2023a, sci prog a3, ex 5
clear; clc;

function [x, w] = quadratureGauss(a, b, n)

    [x01, w01] = gaussLegendre(n);
    x = ((b - a)/2) * x01 + (a + b)/2;
    w = ((b - a)/2) * w01;

end


%% Gauss Legendre Quadrature Weights
% calculate weights and nodes of Gauss Legendre quadrature on [-1,1] with n nodes
% code 'borrowed' from "numerical recipes" (originally in C)
function [x, w] = gaussLegendre(n)
    % Gauss–Legendre Quadrature on [-1,1] with n nodes
    x = zeros(n,1);
    w = zeros(n,1);

    % The roots of the Legendre polynomial P_n are symmetrical in [-1,1].
    % We'll find the roots in [0,1] and mirror them.
    for i = 1:(n+1)/2
        % First approximation for i-th root
        z = cos(pi*(i - 0.25)/(n + 0.5));

        % Newton iteration
        zOld = inf;
        while abs(z - zOld) > 1e-12
            p1 = 1;   % P_0(z)
            p2 = 0;   % P_{-1}(z) dummy

            % Recurrence to compute P_n(z)
            for j = 1:n
                p3 = p2;
                p2 = p1;
                p1 = ((2*j - 1)*z*p2 - (j - 1)*p3)/j;  
            end

            % Derivative of P_n(z) using P_{n-1}(z) = p2
            pp = n*(z*p1 - p2)/(z^2 - 1);

            zOld = z;
            z = z - p1/pp;
        end

        % Set root and its symmetric counterpart
        x(i)         = -z;
        x(n + 1 - i) =  z;

        % Calculate weight (same for both symmetrical roots)
        wVal = 2 / ((1 - z^2)*pp^2);
        w(i)            = wVal;
        w(n + 1 - i)    = wVal;
    end

    % Return row vectors if you prefer
    x = x(:).'; 
    w = w(:).';
end