% kundyz muktar, matlab r2023a, sci prog a3, ex 8
clear; clc; close all;

t = linspace(0,1,100);
a = @(x) ones(size(x));
f = @(x) ones(size(x));

u = p1FEM(a, f, t);
g = 0.5*(1-t).*t;
figure; 
plot(t,u,'.',t,g); 
legend('FEM','Exact');
fprintf('Max error = %.3e\n',max(abs(u-g)));
disp(u)


function u = p1FEM(a, f, t)
%% p1-FEM
% Function to solve p1-FEM in 1D for general second order elliptic PDEs of
% the form
%   -(a(x) u'(x))' + b(x) u'(x) + c(x) u(x) = f(x)
% Input parameters:
%   a   - diffusion coefficient (function handle)
%   f   - right-hand side function (function handle)
%   t   - vector of time steps (triangulation)
% Output:
%   x   - coefficient vector of finite element solution




ndof = numel(t); %???
nIntervals = ndof - 1;


Alocal = zeros(2, 2, nIntervals);
Llocal = zeros(2, nIntervals);




nQuadPoints = 4;

for elem = 1:nIntervals

    xL = t(elem);
    xR = t(elem+1);
    h = xR - xL; %???

    derivatives = [-1/h; 1/h]; %???
    [quadpoints, quadweights] = quadratureGauss(xL, xR, nQuadPoints); %??? we need qg for each el


    % Evaluate local basis functions in local quadrature points
    basisEvaluations = [t(elem+1)-quadpoints; -t(elem)+quadpoints]./h; % basically basis phi for each element
   

    

    % Evaluate integrand of A in local quadrature points
    aEval = a(quadpoints);  
    AEval = zeros(2, 2, nQuadPoints); 

    for j = 1:2
        for k = 1:2
            % same logic: aEval .* derivatives(j)*derivatives(k)
            AEval(j, k, :) = aEval(:) .* (derivatives(j) .* derivatives(k)); %??? by formula
            
        end
    end

    % Compute local stiffness matrix using Gauss quadrature rule
    for k = 1:nQuadPoints %??? 1. should be same as for local rhs 2. stiffness matrix is used in the integral
        Alocal(:, :, elem) = Alocal(:, :, elem) + ...
            quadweights(k) * AEval(:, :, k);

    end

    % Compute local right-hand side
    fEval = f(quadpoints);
    for k = 1:nQuadPoints
        % Llocal is prev defined as a 2:99 matrix of zeros. we fill first
        % row of 1,99 by integrals (as given in the formula)
        Llocal(1, elem) = Llocal(1, elem) + quadweights(k) * fEval(k) * basisEvaluations(1, k);
        

        Llocal(2, elem) = Llocal(2, elem) + quadweights(k) * fEval(k) * basisEvaluations(2, k);
    end

end

% Assemble global stiffness matrix A and load vector L
A = zeros(ndof, ndof);
L = zeros(ndof, 1);

for elem = 1:nIntervals
    % Determine local dofs [1,2],[2,3]..
    localDofs = [elem; elem + 1];
    % Store local stiffness matrix
    % TO DO: Modify the following code line
    %adds the 2×2 local stiffness matrix from the current element into the global matrix 
    % at the rows and columns corresponding to the two nodes of that element.
    A(localDofs, localDofs) = A(localDofs, localDofs) + Alocal(:, :, elem);
    % Store local right-hand side
    L(localDofs) = L(localDofs) + Llocal(:, elem);
end

% Solve linear system of equations
dofs = 2:(ndof - 1);
x = zeros(ndof, 1);
% TO DO: Modify the following code line to compute the FEM solution x
x(dofs) = A(dofs, dofs) \ L(dofs);

u = x;


end


%%
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