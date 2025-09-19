%% explicitRungeKuttaSystem.m
clc; clear; close all;
% ____________________params________________________________
t0 = 0;
tEnd = 1;
N = 10;  % number of subintervals
t = linspace(t0, tEnd, N+1);
% ____________________problem definitions________________________________

d = 5; % dim
f1 = @(tt) zeros(d,1); % d dim vect
f2 = @(tt) -eye(d);% dxd matrix
y0 = ones(d,1);  %ic

a = [0    0    0    0
     1/2  0    0    0
     0   1/2   0    0
     0    0    1    0 ];
b = [1/6  1/3  1/3  1/6];
c = [0; 1/2; 1/2; 1];

% ______________________solving______________________________
y_rk4 = explicitRungeKuttaSystem(a, b, c, f1, f2, y0, t);
y_exact = exp(-t);

% ____________________displaying sols________________________________
disp('time grid t:');
disp(t);
disp('y_rk4 (each column is the solution at time steps):');
disp(y_rk4);
disp('y_exact:');
disp(y_exact);

% _______________________plotting_____________________________
figure; hold on;
plot(t, y_exact, 'm-', 'displayname','exact');
plot(t, y_rk4, 'ro--', 'displayname','rk4');
xlabel('t');
ylabel('y(t)');
title('ivp: y''(t) = -y(t), y(0)=1, extended to dimension d=1');
legend('location','best');
grid on;
hold off;

%% local function definition
function Y = explicitRungeKuttaSystem(a, b, c, f1, f2, y0, t)

L = length(t) - 1; % number of steps
d = length(y0); % dimension
Y = zeros(d, L+1); % store the solution
Y(:,1) = y0; % initial condition
s = length(b); % number of stages

for n = 1:L
    h = t(n+1) - t(n); % step size
    K = zeros(d, s); % each stage is a d dimvector
    for j = 1:s
        % partial sum of previous stages
        sumPrev = zeros(d,1);
        for m = 1:(j-1)
            sumPrev = sumPrev + a(j,m)*K(:,m);
        end
        tj = t(n) + c(j)*h;   % stage time
        yj = Y(:,n) + h*sumPrev;  % stage state
        % evaluate f(tj, yj) = f1(tj) + f2(tj)*yj
        K(:,j) = f1(tj) + f2(tj)*yj;
    end
    % combine stages
    sumStages = zeros(d,1);
    for j = 1:s
        sumStages = sumStages + b(j)*K(:,j);
    end
    % final update
    Y(:,n+1) = Y(:,n) + h*sumStages;
end
end
