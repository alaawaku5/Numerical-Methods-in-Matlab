Scientific Programming, Worksheet 4, Ex. 3
Written by: Kundyz Muktar, Course: Scientific Programming (2025), MATLAB version: R2023a




Taken from Worksheet 3, Ex 3:
f = @(x)exp(-x);
f_exact = 23.07470459693395647;



err_vals = zeros(1,20);
for i=1:20

    [x_vals, w] = quadrature(-pi, exp(1),i);
    [A, b_vals] = build_A_b(i-1,i,x_vals);

    Q = sum(w .* f(x_vals)');
    err_vals(i) = abs(Q - f_exact);

end


figure;
semilogy(1:20, err_vals, '-o');
xlabel('# of quadrature nodes (n)');
ylabel('absolute error');



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
