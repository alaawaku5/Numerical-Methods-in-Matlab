% kundyz muktar, matlab r2023a, sci prog a3, ex 7, part 2
clear; clc;

f = @(x) 1./(1+(x.^2));
f_exact = 2.74680153389003172;




err_vals = zeros(1,50);
for i=1:50

    [x_vals, w] = quadrature(-5, 5,i);
    [A, b_vals] = build_A_b(i-1,i,x_vals);

    Q = sum(w .* f(x_vals)');
    err_vals(i) = abs(Q - f_exact);

end
disp(Q);

figure;
semilogy(1:50, err_vals, '-o');
xlabel('# of quadrature nodes (n)');
ylabel('absolute error');




%% from ex 3 
function [x, w] = quadrature(a, b, n)
    if n == 1
        x = a;  % if only one node, choose the left endpoint (could also choose the midpoint)
        w = (b-a);
        return
    end

    h = (b - a) / (n - 1);  % quadrature step size
    x = a:h:b;
    

    x01 = linspace(0, 1, n); 
    w01 = quadratureWeights(x01);
    w = (b - a) * w01;
end

%% from ex 1 or 2
function w = quadratureWeights(x)
    n = numel(x);
    A = zeros(n,n);
    b = zeros(n,1);

    for row = 0:n-1
        A(row+1,:) = x.^row;
        b(row+1)   = 1/(row+1);
    end

    w = A\b;
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