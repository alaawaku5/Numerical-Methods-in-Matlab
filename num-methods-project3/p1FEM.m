function x = p1FEM(a, f, t)
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

% TO DO: Modify the following code line
ndof = %???
nIntervals = ndof - 1;

% Initialise arrays
Alocal = zeros(2, 2, nIntervals);
Llocal = zeros(2, nIntervals);

% Number of employed quadrature points
nQuadPoints = 4;

% Compute local contributions to stiffness matrix S and right-hand side F
for elem = 1:nIntervals
    % Compute length of interval
    % TO DO: Modify the following code line
    h = %???

    % Compute derivatives of basis functions as column vector
    % TO DO: Modify the following code line
    derivatives = %???

    % Compute the Gauss quadrature on the current interval
    % with the function quadratureGauss(a, b, nQuadPoints), where a and b
    % are the left and right endpoints, respectively, and the number of
    % quadrature points nQuadPoints
    % TO DO: Modify the following code line
    [quadpoints, quadweights] = quadratureGauss(%???, %???, nQuadPoints);

    % Evaluate local basis functions in local quadrature points
    basisEvaluations = [t(elem+1)-quadpoints; -t(elem)+quadpoints]./h;

    % Evaluate integrand of A in local quadrature points
    aEval = a(quadpoints);
    AEval = zeros(2, 2, nQuadPoints);
    for j = 1:2
        for k = 1:2
            % TO DO: Modify the following code line to evaluate
            % the integrand of A(j,k) in the quadrature points
            % as in equation (6)
            AEval(j, k, :) = %???
        end
    end

    % Compute local stiffness matrix using Gauss quadrature rule
    % TO DO: Find out the maximal index for the k-loop
    for k = 1:%???
        Alocal(:, :, elem) = Alocal(:, :, elem) + ...
            quadweights(k) * AEval(:, :, k);
    end


    % Compute local right-hand side
    fEval = f(quadpoints);
    for k = 1:nQuadPoints
        Llocal(1, elem) = Llocal(1, elem) + ...
            quadweights(k) * fEval(k) * basisEvaluations(1, k);
        % TO DO: Modify the following code line
        Llocal(2, elem) = %???
    end
end

% Assemble global stiffness matrix S and load vector F
A = zeros(ndof, ndof);
L = zeros(ndof, 1);
for elem = 1:nIntervals
    % Determine local dofs
    localDofs = [elem; elem + 1];
    % Store local stiffness matrix
    % TO DO: Modify the following code line
    A(localDofs, localDofs) = %???;
    % Store local right-hand side
    L(localDofs) = L(localDofs) + Llocal(:, elem);
end

% Solve linear system of equations
dofs = 2:(ndof-1);
x = zeros(ndof, 1);
% TO DO: Modify the following code line to compute the FEM solution x
x(dofs) = %???

end
