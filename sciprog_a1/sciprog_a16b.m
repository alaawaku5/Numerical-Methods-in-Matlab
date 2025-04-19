% phaseFieldPlot_2DSystemScript.m
% plotting the phase field (y2 vs y1) for the 2d system y1'(t)=y2(t), y2'(t)=-y1(t)

clc; clear; close all;

h = 0.1; % step size
t = 0:h:(2*pi); % time grid

d = 2; % dimension of system
f1 = @(tt) zeros(d,1); % zero forcing
f2 = @(tt) [0 1; -1 0]; % matrix for system

y0 = [1; 0]; % initial condition

% solve with each method
y_explicit = explicit_Euler(t, f1, f2, y0);
y_implicit = implicit_Euler(t, f1, f2, y0);
y_midpoint = implicit_Midpoint(t, f1, f2, y0);

% exact solution on same grid
y_exact = zeros(d, length(t));
for n = 1:length(t)
    y_exact(:,n) = [cos(t(n)); sin(t(n))];
end

% phase field plot: x-axis = y2(t), y-axis = y1(t)
figure; hold on;
plot(y_exact(2,:), y_exact(1,:), 'o-', 'displayname','exact'); % exact sol
plot(y_explicit(2,:), y_explicit(1,:), 'm', 'displayname','explicit euler'); % explicit
plot(y_implicit(2,:), y_implicit(1,:), 'g', 'displayname','implicit euler'); % implicit
plot(y_midpoint(2,:), y_midpoint(1,:), 'r', 'displayname','implicit midpoint'); % midpoint

xlabel('y_2(t)'); % x axis
ylabel('y_1(t)'); % y axis
title('phase field plot with h=0.1'); % short title
legend('location','best'); % show legend
axis equal; % keep circle shape
grid on;
hold off;

%% local function definitions

function y = explicit_Euler(t, f1, f2, y0)
d = length(y0); % dimension
L = length(t) - 1; % steps
y = zeros(d, L+1); % store sol
y(:,1) = y0; % init
for n = 1:L
    h = t(n+1) - t(n); % step
    A = f1(t(n)); % f1
    B = f2(t(n)); % f2
    y(:,n+1) = y(:,n) + h*(A + B*y(:,n)); % update
end
end

function y = implicit_Euler(t, f1, f2, y0)
d = length(y0); % dimension
L = length(t) - 1; % steps
y = zeros(d, L+1); % store sol
y(:,1) = y0; % init
for n = 1:L
    h = t(n+1) - t(n); % step
    A = f1(t(n+1)); % f1
    B = f2(t(n+1)); % f2
    I = eye(d); % identity
    y(:,n+1) = (I - h*B)\(y(:,n) + h*A); % solve system
end
end

function y = implicit_Midpoint(t, f1, f2, y0)
d = length(y0); % dimension
L = length(t) - 1; % steps
y = zeros(d, L+1); % store sol
y(:,1) = y0; % init
for n = 1:L
    h = t(n+1) - t(n); % step
    tm = (t(n) + t(n+1))/2; % midpoint
    A = f1(tm); % f1
    B = f2(tm); % f2
    I = eye(d); % identity
    y(:,n+1) = (I - (h/2)*B)\( y(:,n) + h*A + (h/2)*B*y(:,n) ); % solve
end
end
