%% ex 9
%   (a) Compute the error and compare the convergence rate with the results from Problem 6 (a).
% What convergence rate do you observe?
%   (b) Plot the solution into the phase field plot from Problem 6 (b).

clc; clear; close all;


d = 2;  % dim of the system

if d == 2
    f1 = @(tt) zeros(d,1);    % 2x1 zero vector
    f2 = @(tt) [0 1; -1 0];   % 2x2
    y0 = [1; 0];  % initial condition
    exact_sol = @(tt) [cos(tt); sin(tt)]; % for error check
else
    % 1d ex
    f1 = @(tt) zeros(d,1);   
    f2 = @(tt) -eye(d);    
    y0 = ones(d,1);
    exact_sol = @(tt) exp(-tt)*ones(d,1);
end

% butcher tableau for classic 4-stage rk4
a = [0    0    0    0
     1/2  0    0    0
     0   1/2   0    0
     0    0    1    0 ];
b = [1/6  1/3  1/3  1/6];
c = [0; 1/2; 1/2; 1];

% PART A
% we'll vary step sizes h = 2^-i for i=1..L
L = 10;  % up to i=5 => h=1/2,1/4,1/8,1/16,1/32
error_vals = zeros(1,L);
h_vals = zeros(1,L);

for i = 1:L
    h = 2^(-i);
    h_vals(i) = h;
    tEnd = 1;    % final time
    t = 0:h:tEnd;
    % lets solve
    yrk = explicitRungeKuttaSystem(a, b, c, f1, f2, y0, t);

    % build exact sol on same grid
    y_exactact = zeros(d, length(t));
    for n=1:length(t)
        y_exact(:,n) = exact_sol(t(n));
    end

  % ******* ERROR *******
    % measure error in inf norm over all time steps, max error in each
    % vector in each step, and the max error across all steps
    diff = yrk - y_exact;
    stepNorms = vecnorm(diff, Inf, 1);
    error_vals(i) = max(stepNorms);
end

% compute approximate order
orderVals = zeros(1, L-1);
for j = 1:L-1
    orderVals(j) = log(error_vals(j)/error_vals(j+1)) / log(2);
end

% show results
disp('Step sizes h_vals:');
disp(h_vals);
disp('Error values (inf norm):');
disp(error_vals);
disp('Estimated order of convergence (pairwise):');
disp(orderVals);

% log-log plot of error
figure;
loglog(h_vals, error_vals, 'ro--'); 
xlabel('step size h');
ylabel('error (inf norm)');
title('Error vs. Step size (log-log)');
grid on;

% Part B 
if d == 2
    % pick a single step size to see the shape, e.g. h=0.1
    h_size = 0.1;
    t_vals = 0:h_size:2*pi; % go to 2*pi to see a full cycle
    y_approx = explicitRungeKuttaSystem(a, b, c, f1, f2, y0, t_vals);

    % build exact on same grid
    y_exactact = zeros(2, length(t_vals));
    for n=1:length(t_vals)
        y_exactact(:,n) = exact_sol(t_vals(n)); % exact sol at each step
    end

    figure; hold on;
    % problem 6 used y2 vs. y1, so let's do that:
    plot(y_exactact(1,:), y_exactact(2,:), 'k-', 'displayname','exact');
    plot(y_approx(1,:), y_approx(2,:), 'ro--', 'displayname','rk4');
    xlabel('y_1(t)');
    ylabel('y_2(t)');
    title('phase field plot (d=2) with h=0.1');
    axis equal; % keep shape
    legend('location','best');
    grid on;
    hold off;
end
%% more rapid convergence in this
%% local function
function y = explicitRungeKuttaSystem(a, b, c, f1, f2, y0, t)
L = length(t) - 1; % # steps
d = length(y0); % # system dum
y = zeros(d, L+1); % sol in d by L+1 matrix
y(:,1) = y0; % ic
s = length(b); % number of stages

for n = 1:L
    h = t(n+1) - t(n); % size of the sep
    K = zeros(d, s); % storing STAGE derivative, for each s we have d sized vector
    %*********stage***********
    for j = 1:s
        sumPrev = zeros(d,1);
        for m = 1:j-1
            sumPrev = sumPrev + a(j,m)*K(:,m); % a here is from the given
        end
        tj = t(n) + c(j)*h; % time for jth stage 
        yj = y(:,n) + h*sumPrev; %
        K(:,j) = f1(tj) + f2(tj)*yj;
    end
    %********************
    sumStages = zeros(d,1);
    for j = 1:s
        sumStages = sumStages + b(j)*K(:,j);
    end
    y(:,n+1) = y(:,n) + h*sumStages;
end
end
