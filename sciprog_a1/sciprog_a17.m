%% Explicit Runge Kutta
% without a doubt, runge kutta has a higher accuracy than the euler methods
clc; clear; close all;

t0 = 0;
tEnd = 1;
N = 10; % # subintervals
t = linspace(t0, tEnd, N+1); % equidistant grid

% for y'(t)=-y(t), we have f1(t)=0, f2(t)=-1
f1 = @(tt) 0;       
f2 = @(tt) -1;       
y0 = 1; % initial condition

a = [0    0    0    0
     1/2  0    0    0
     0   1/2   0    0
     0    0    1    0 ];
b = [1/6  1/3  1/3  1/6];
c = [0; 1/2; 1/2; 1];

% solve system via rk4 and exact
y_rk4 = explicitRungeKutta(a, b, c, f1, f2, y0, t);
y_exact = exp(-t);

% --- display results
disp('time grid t:');
disp(t);
disp('y rk4:');
disp(y_rk4);
disp('y_exact:');
disp(y_exact);

% --- plot
figure; hold on;
plot(t, y_exact, 'm-', 'DisplayName','exact');         % exact
plot(t, y_rk4, 'ro--', 'DisplayName','rk4 numeric');   % numeric
xlabel('t');
ylabel('y(t)');
title('IVP: y''(t)=-y(t), y(0)=1');
legend('Location','best');
grid on;
hold off;

%% local function definition
function y = explicitRungeKutta(a, b, c, f1, f2, y0, t)
L = length(t) - 1; % # steps
y = zeros(1, L+1);  % store solution for future retrieval
y(1) = y0; % initial condition
s = length(b); % # stages

for n = 1:L
    h = t(n+1)-t(n);   % step size
    K = zeros(1, s);     % each k_j is scalar in this problem
    for j = 1:s
        % partial sum for the previous stages
        sumPrev = 0;
        for m = 1:(j-1)
            sumPrev = sumPrev + a(j,m)*K(m);
        end
        tj = t(n) + c(j)*h; % time for stage j
        yj = y(n) + h*sumPrev;  % approximate y for stage j
        K(j) = f1(tj) + f2(tj)*yj; % evaluate derivative: k_j = f1(tj) + f2(tj)*yj
    end
    % combine stages
    sumStages = 0;
    for j = 1:s
        sumStages = sumStages + b(j)*K(j);
    end
    % final update
    y(n+1) = y(n) + h*sumStages;
end
end
