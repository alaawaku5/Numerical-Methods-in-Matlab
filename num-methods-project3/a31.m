% kundyz muktar, matlab r2023a, sci prog a3, ex 1
clear; clc;
N = 5;    % quadrature points
k = 4;    % exactness up to x^k

h = 1/N;
x_vals = (0:h:1-h);         % 1×N row vector



[A,b_vals] = build_A_b(k,N,x_vals);

%%%%%%%%%%%%%%%%%%%(Part A)%%%%%%%%%%%%%%%%%%%
w = A \ b_vals;
% isvector(w)
disp(w)
disp((10)*w)


%%%%%%%%%%%%%%%%%%%(Part B)%%%%%%%%%%%%%%%%%%%
% n point quadrature is exact up to order n-1, 
% so we can take the degrees of pols up to n, and have a sqr matrix



function w = quadratureWeights(x)
    n = numel(x); 
    A = zeros(n,n);
    b = zeros(n, 1);
    
    % fill up A by powers of x
    for row =0:n-1
        A(row+1,:) = x.^row;
        b(row+1) = 1/(row+1);
    end
    w = A \ b;

end

%%
% used for part A
function [A,b_vals] = build_A_b(k,N,x_vals)

    A = ones(k+1, N);
    b_vals = zeros(k+1, 1);
    for row = 0:k
        A(row+1, :) = x_vals.^row;
        b_vals(row+1) = 1/(row+1);
    end
end