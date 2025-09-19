% kundyz muktar, matlab r2023a, sci prog a3, ex 9-10
clear; clc; close all;

a = @(x) 1;                   
c = @(x) -1;                  
f = @(x) ones(size(x));       

bVals = [-10, -5, 0, 5, 10];

n = 100;                   
t = linspace(0, 1, n+1);   

% solve for each b and plot --
figure('Name','Diffusion-Advection-Reaction');
hold on;
colors = lines(length(bVals));  % for nicer plotting
%  b is the advection term, in some problems the advection term 
% represents a physical transport mechanism
for i = 1:length(bVals)
    bConst = bVals(i);
    b = @(x) bConst;

    xSol = p1FEM(a, b, c, f, t);

    plot(t, xSol, 'Color', colors(i,:),'DisplayName', sprintf('b = %d', bConst));
end

xlabel('x');
ylabel('u(x)');
title('Solutions to -(u'') + b*u'' - u = 1 for various b');
legend('Location','best');
grid on;
hold off;


function u = p1FEM(a, b, c, f, t)
%% p1-FEM for the PDE
%   -(a(x) u'(x))' + b(x) u'(x) + c(x) u(x) = f(x),  on [0,1]
% with Dirichlet BCs  u(0)=0, u(1)=0.
%
% INPUT:
%   a   - function handle for diffusion coefficient  a(x)
%   b   - function handle for advection coefficient b(x)
%   c   - function handle for reaction coefficient  c(x)
%   f   - function handle for right-hand side       f(x)
%   t   - row vector of node coordinates (mesh), length(t)=N+1
%
% OUTPUT:
%   u   - finite element solution at each mesh node

    ndof = numel(t);        % total number of nodes
    nIntervals = ndof - 1;  % number of elements

    % Local containers for each element
    Alocal = zeros(2, 2, nIntervals);  % local 2x2 stiffness for each elem
    Llocal = zeros(2, nIntervals);     % local 2x1 load for each elem

    nQuadPoints = 4;  % e.g. 4-point Gauss rule

    for elem = 1:nIntervals

        % Element interval [xL, xR] and length
        xL = t(elem);
        xR = t(elem+1);
        h  = xR - xL;

        % Derivatives of the 2 local shape functions (piecewise linear)
        %   phi1'(x) = -1/h,   phi2'(x) = +1/h
        derivatives = [-1/h; 1/h];

        % Quadrature points & weights on [xL, xR]
        [quadpoints, quadweights] = quadratureGauss(xL, xR, nQuadPoints);

        % Evaluate shape functions at the quadrature points:
        %   phi1(x) = (xR - x)/h
        %   phi2(x) = (x - xL)/h
        basisEvaluations = [ (xR - quadpoints);
                             (quadpoints - xL) ] ./ h;   % 2 x nQuadPoints
        basisEvaluations
        % Evaluate PDE coefficients at quadpoints
        aEval = a(quadpoints);    % diffusion
        bEval = b(quadpoints);    % advection
        cEval = c(quadpoints);    % reaction

        % Build integrand for local stiffness:
        % A = ∫ [ a(x)*phi_i'(x)*phi_j'(x)
        %       + b(x)*phi_i'(x)*phi_j(x)
        %       + c(x)*phi_i(x)*phi_j(x) ] dx
        AEval = zeros(2, 2, nQuadPoints);
        for i = 1:2
            for j = 1:2
                % a(x)*dphi_i*dphi_j
                aPart = aEval(:) .* derivatives(i) .* derivatives(j);

                % b(x)*dphi_i*phi_j
                bPart = bEval(:) .* derivatives(i) .* basisEvaluations(j,:);

                % c(x)*phi_i*phi_j
                cPart = cEval(:) .* basisEvaluations(i,:) .* basisEvaluations(j,:);

                % Sum them up for each quadrature point
                AEval(i, j, :) = aPart + bPart + cPart;
            end
        end

        % Integrate AEval over the element via Gauss quadrature
        for q = 1:nQuadPoints
            Alocal(:,:,elem) = Alocal(:,:,elem) + ...
                quadweights(q) * AEval(:,:,q);
        end

        % Build local load vector: L = ∫ f(x)*phi_i(x) dx
        fEval = f(quadpoints);
        for q = 1:nQuadPoints
            % contribution to phi1
            Llocal(1, elem) = Llocal(1, elem) + ...
                quadweights(q) * fEval(q) * basisEvaluations(1, q);

            % contribution to phi2
            Llocal(2, elem) = Llocal(2, elem) + ...
                quadweights(q) * fEval(q) * basisEvaluations(2, q);
        end

    end  % end element loop

    %--- Global assembly of A and L ---
    A = zeros(ndof, ndof);
    L = zeros(ndof, 1);

    for elem = 1:nIntervals
        localDofs = [elem; elem+1];
        % Insert 2x2 local matrix into global matrix
        A(localDofs, localDofs) = A(localDofs, localDofs) + Alocal(:,:,elem);
        % Insert 2x1 local load into global load
        L(localDofs) = L(localDofs) + Llocal(:,elem);
    end

    %--- Impose Dirichlet boundary conditions u(0)=0, u(1)=0 ---
    dofs = 2:(ndof-1);
    x = zeros(ndof, 1);

    % Solve reduced system
    x(dofs) = A(dofs, dofs)\L(dofs);

    % Return the full solution vector
    u = x;
end



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
