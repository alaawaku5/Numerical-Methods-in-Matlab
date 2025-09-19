function [x,w]=gaussLegendre(n)
%% Gauss Legendre Quadrature Weights
% calculate weights and nodes of Gauss Legendre quadrature on [-1,1] with n nodes
% code 'borrowed' from "numerical recipes" (originally in C)
%
% Input parameters:
%   n   - number of nodes
%
% Output:
%   x   - quadrature nodes
%   w   - corresponding quadrature weights
%
% This function is part of the course
%
%   "Scientific programming for interdisciplinary mathematics"
%   Summer 2025
%
%   Dr. Markus Faustmann
%   Institute of Analysis and Scientific Computing
%   TU Wien

    x = zeros(n,1);
    w = zeros(n,1);

    % The roots of the Legendre polynomial P_n are symmetrical in [-1,1],
    %   i.e., we only need to find half of them
    for i=1:(n+1)/2
        % first approximation of i-th root
        z = cos(pi*(i - 0.25)/(n + 0.5));

        % Newton's method
        zOld = inf;
        while abs(z - zOld) > 1e-12
            % calculate p1 = P_n(z) by the recursion formula
            %   (k+1)*P_{k+1}(x) = (2k+1)*x*P_k(x) - k*P_{k-1}(x)
            p1 = 1;
            p2 = 0;
            for j=1:n
                p3 = p2;
                p2 = p1;
                p1 = ((2*j - 1)*z*p2 - (j - 1)*p3)/j;
            end

            % calculate pp = P'_n(z) using p2 = P_{n-1}(z)
            pp = n*(z*p1 - p2)/(z*z - 1);

            % update z
            zOld = z;
            z = z - p1/pp;
        end

        % set root and it's symmetric counterpart
        x(i) = -z;
        x(n+1-i) = z;

        % calculate weight and set symmetric counterpart
        w(i) = 2 / ((1 - z^2)*pp^2);
        w(n+1-i) = w(i);
    end
    x = x';
end
