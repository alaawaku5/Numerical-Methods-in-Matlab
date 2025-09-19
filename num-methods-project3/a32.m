% kundyz muktar, matlab r2023a, sci prog a3, ex 2
%CHECK WHY LOGLOG PLOT IS NOT STRAIGHT LINE
clear; clc;

%%%%%%%%%%%%%%%%%%%(Part 1)%%%%%%%%%%%%%%%%%%%
n = 20;

%%%%%%%%%%%%%%%%%%%(Part 2)%%%%%%%%%%%%%%%%%%%
f1 = @(x) exp(x);
f2 = @(x) x.^(3/2);

err_vals1 = zeros(1,n);
err_vals2 = zeros(1,n);

for i = 1:n
    x = (0:i)/i;                       % <-- correct node definition
    weights = quadratureWeights(x);    % <-- recompute weights each i
    
    Q1 = sum(weights .* f1(x)'); 
    Q2 = sum(weights .* f2(x)');


    err_vals1(i) = abs(Q1 - (exp(1)-1));
    err_vals2(i) = abs(Q2 - 2/5);
end

disp(err_vals2)

nVec = 1:n;

% f(x)=exp(x), exponential decay?
figure
semilogy(nVec, err_vals1, '-o')
xlabel('number of nodes n')
ylabel('abs error')
title('semilogy, error_{n} for f(x)=exp(x)')

% f(x)=x^{3/2} expect algebraic decay?
figure
loglog(nVec, err_vals2, '-o')
xlabel('Number of nodes n')
ylabel('Absolute error')
title('loglog, error_{n} for f(x)=x^{3/2}')


%% from ex 1
function w = quadratureWeights(x)
    n = numel(x);
    A = ones(n,n);
    b = zeros(n,1);

    for row = 0:n-1
        A(row+1,:) = x.^row;
        b(row+1)   = 1/(row+1);
    end

    w = A\b;
end
