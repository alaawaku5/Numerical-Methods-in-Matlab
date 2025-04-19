% kundyz muktar, matlab r2023a, sci prog a3, ex 6
clear; clc;



f = @(x)exp(-x);
f_exact = 23.07470459693395647;


err_vals = zeros(1,20);
for i=1:20

    [x_vals, w] = quadratureGauss(-pi, exp(1),i); % only this part is diff from ex 3
    [A, b_vals] = build_A_b(i-1,i,x_vals);

    Q = sum(w .* f(x_vals));
    err_vals(i) = abs(Q - f_exact);

end
disp(Q);

figure;
semilogy(1:20, err_vals, '-o');
xlabel('# of quadrature nodes (n)');
ylabel('absolute error');







%% from ex 5
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

%% from ex 1
function [A, b_vals] = build_A_b(k, N, x_vals)
    % Here k is the maximum polynomial degree (so there are k+1 moments)
    A = zeros(k+1, N);
    b_vals = zeros(k+1, 1);
    for row = 0:k
        A(row+1, :) = x_vals.^row;
        b_vals(row+1) = 1/(row+1);
    end
end