% kundyz muktar, matlab r2023a, sci prog a3, ex 8
clear; clc; close all;

% Script code for testing (EXERCISE 8 PART, LEFT INTACT)
a = @(x) 1;
f = @(x) ones(size(x));
n = 100;
t = linspace(0, 1, n+1);
xSol = p1FEM(a, f, t);

g = @(x) 0.5 .* x .* (1 - x);
figure;
plot(t, xSol, 'ro-', 'LineWidth', 1.2, 'DisplayName', 'FEM solution'); 
hold on;
plot(t, g(t), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Exact solution');
grid on;
xlabel('x');
ylabel('u(x)');
title('Comparison: p1FEM solution vs. exact');
legend('Location','best');

%
%%
function x = p1FEM(a, f, t)
    ndof = length(t); %???
    nIntervals = ndof - 1;
    
    Alocal = zeros(2, 2, nIntervals);
    Llocal = zeros(2, nIntervals);
    
    nQuadPoints = 4;
    
    for elem = 1:nIntervals
    
        % Element length
        h = t(elem+1) - t(elem); %???
    
        % Derivatives of the local shape functions
        derivatives = [-1/h; 1/h]; %???
    
        % Quadrature points and weights on [t(elem), t(elem+1)]
        [quadpoints, quadweights] = quadratureGauss(t(elem), t(elem+1), nQuadPoints); %???
    
        % Evaluate the basis functions at each quadrature point
        %  phi1(x) = (t(elem+1)-x)/h,   phi2(x) = (x - t(elem))/h
        basisEvaluations = [t(elem+1)-quadpoints; -t(elem)+quadpoints]./h;
    
        % Evaluate coefficient a(x) at quadrature points
        aEval = a(quadpoints);
    
        % Precompute local integrand for stiffness
        AEval = zeros(2, 2, nQuadPoints);
        for j = 1:2
            for k = 1:2
                AEval(j, k, :) = aEval(:) .* derivatives(j) .* derivatives(k); %???
            end
        end
    
        % Integrate to form Alocal
        for k = 1:nQuadPoints %???
            Alocal(:, :, elem) = Alocal(:, :, elem) + quadweights(k) * AEval(:, :, k);
        end
    
        % Integrate to form Llocal
        fEval = f(quadpoints);
        for k = 1:nQuadPoints
            Llocal(1, elem) = Llocal(1, elem) + quadweights(k) * fEval(k) * basisEvaluations(1, k);
            Llocal(2, elem) = Llocal(2, elem) + quadweights(k) * fEval(k) * basisEvaluations(2, k); %???
        end
    end
    
    % Global assembly
    A = zeros(ndof, ndof);
    L = zeros(ndof, 1);
    for elem = 1:nIntervals
        localDofs = [elem; elem + 1];
        A(localDofs, localDofs) = A(localDofs, localDofs) + Alocal(:, :, elem); %???
        L(localDofs) = L(localDofs) + Llocal(:, elem);
    end
    
    % Solve for the interior dofs (Dirichlet BCs at t(1) and t(end))
    dofs = 2:(ndof-1);
    x = zeros(ndof, 1);
    x(dofs) = A(dofs, dofs)\L(dofs); %???
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% additional code for problem 11 (tracking computational time and sparse usage)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% small letters for comments:
% we will measure cpu times for n=2^k and compare full vs. sparse approach

% define the range of n=2^k
kVals = 4:10; 
nVals = 2.^kVals;

timeFull   = zeros(size(nVals));
timeSparse = zeros(size(nVals));

% loop over n=2^k
for i = 1:length(nVals)
    nTest = nVals(i);
    tTest = linspace(0, 1, nTest+1);

    % measure time for full approach
    tic
    xFull = p1FEM_Full(a, f, tTest);
    timeFull(i) = toc;

    % measure time for sparse approach
    tic
    xSparse = p1FEM_Sparse(a, f, tTest);
    timeSparse(i) = toc;
end

% plot results on a loglog scale
figure;
loglog(nVals, timeFull, 'ro-', 'LineWidth', 1.5, 'DisplayName','full');
hold on;
loglog(nVals, timeSparse, 'bs-', 'LineWidth', 1.5, 'DisplayName','sparse');
grid on;
xlabel('n (number of elements)');
ylabel('time (sec)');
title('comparison: full vs. sparse');
legend('Location','best');


% small letters for comments:
% p1fem_full is the same as p1fem but we keep matrix in full format
function x = p1FEM_Full(a, f, t)
    ndof = length(t);
    nIntervals = ndof - 1;
    
    Alocal = zeros(2,2,nIntervals);
    Llocal = zeros(2,nIntervals);

    nQuadPoints = 4;

    for elem = 1:nIntervals
        xL = t(elem);
        xR = t(elem+1);
        h  = xR - xL;

        derivatives = [-1/h; 1/h];
        [quadpoints, quadweights] = quadratureGauss(xL, xR, nQuadPoints);
        basisEvals = [xR - quadpoints; quadpoints - xL]./h;
        
        aEval = a(quadpoints);
        
        AEval = zeros(2,2,nQuadPoints);
        for iLoc = 1:2
            for jLoc = 1:2
                AEval(iLoc, jLoc, :) = aEval .* derivatives(iLoc) .* derivatives(jLoc);
            end
        end

        for q = 1:nQuadPoints
            Alocal(:,:,elem) = Alocal(:,:,elem) + quadweights(q)*AEval(:,:,q);
        end

        fEval = f(quadpoints);
        for q = 1:nQuadPoints
            Llocal(1,elem) = Llocal(1,elem) + quadweights(q)*fEval(q)*basisEvals(1,q);
            Llocal(2,elem) = Llocal(2,elem) + quadweights(q)*fEval(q)*basisEvals(2,q);
        end
    end

    A = zeros(ndof, ndof);
    L = zeros(ndof,1);

    for elem = 1:nIntervals
        localDofs = [elem; elem+1];
        A(localDofs, localDofs) = A(localDofs, localDofs) + Alocal(:,:,elem);
        L(localDofs)           = L(localDofs)           + Llocal(:,elem);
    end

    dofs = 2:(ndof-1);
    x = zeros(ndof,1);
    x(dofs) = A(dofs,dofs)\L(dofs);
end


% small letters for comments:
% p1fem_sparse uses a sparse global matrix for the same problem
function x = p1FEM_Sparse(a, f, t)
    ndof = length(t);
    nIntervals = ndof - 1;
    
    Alocal = zeros(2,2,nIntervals);
    Llocal = zeros(2,nIntervals);

    nQuadPoints = 4;

    for elem = 1:nIntervals
        xL = t(elem);
        xR = t(elem+1);
        h  = xR - xL;

        derivatives = [-1/h; 1/h];
        [quadpoints, quadweights] = quadratureGauss(xL, xR, nQuadPoints);
        basisEvals = [xR - quadpoints; quadpoints - xL]./h;
        
        aEval = a(quadpoints);
        
        AEval = zeros(2,2,nQuadPoints);
        for iLoc = 1:2
            for jLoc = 1:2
                AEval(iLoc, jLoc, :) = aEval .* derivatives(iLoc) .* derivatives(jLoc);
            end
        end

        for q = 1:nQuadPoints
            Alocal(:,:,elem) = Alocal(:,:,elem) + quadweights(q)*AEval(:,:,q);
        end

        fEval = f(quadpoints);
        for q = 1:nQuadPoints
            Llocal(1,elem) = Llocal(1,elem) + quadweights(q)*fEval(q)*basisEvals(1,q);
            Llocal(2,elem) = Llocal(2,elem) + quadweights(q)*fEval(q)*basisEvals(2,q);
        end
    end

    % we store triplets for the global matrix
    rowInds = zeros(4*nIntervals,1);
    colInds = zeros(4*nIntervals,1);
    vals    = zeros(4*nIntervals,1);
    L = zeros(ndof,1);

    idx = 1;
    for elem = 1:nIntervals
        localDofs = [elem; elem+1];
        L(localDofs) = L(localDofs) + Llocal(:,elem);

        for iLoc = 1:2
            for jLoc = 1:2
                rowInds(idx) = localDofs(iLoc);
                colInds(idx) = localDofs(jLoc);
                vals(idx)    = Alocal(iLoc,jLoc,elem);
                idx = idx+1;
            end
        end
    end

    A = sparse(rowInds, colInds, vals, ndof, ndof);

    dofs = 2:(ndof-1);
    x = zeros(ndof,1);
    x(dofs) = A(dofs,dofs)\L(dofs);
end
