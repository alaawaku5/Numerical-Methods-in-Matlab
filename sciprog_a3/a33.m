% kundyz muktar, matlab r2023a, sci prog a3, ex 3
clear; clc;

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