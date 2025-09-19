clc; clear; close all;

t0 = 0; tEnd = 1;
N = 80; % # of subintervals
t = linspace(t0, tEnd, N+1);
%--------------------------------------------------
% y'(t)=-y(t), y(0)=1
f1 = @(tt) 0; % f1=0
f2 = @(tt) -1; 
y0 = 1;
y_exact = exp(-t);
%----------------2-stage gauss----------------------------------
a2 = [1/4,    1/4 - sqrt(3)/6;
      1/4 + sqrt(3)/6,   1/4];
b2 = [1/2, 1/2];
c2 = [1/2 - sqrt(3)/6; 1/2 + sqrt(3)/6];
%-------------------3-stage gauss-------------------------------
a3 = [...
   5/36,                 2/9 - sqrt(15)/15,   5/36 - sqrt(15)/30
   5/36 + sqrt(15)/24,   2/9,                 5/36 - sqrt(15)/24
   5/36 + sqrt(15)/30,   2/9 + sqrt(15)/15,   5/36               ];
b3 = [5/18, 4/9, 5/18];
c3 = [1/2 - sqrt(15)/10; 1/2; 1/2 + sqrt(15)/10];
%--------------------------------------------------
gauss2_sol = implicitRK_scalarLinear(a2, b2, c2, t, f1, f2, y0);

gauss3_sol = implicitRK_scalarLinear(a3, b3, c3, t, f1, f2, y0);
%--------------------------------------------------
% display
disp('time grid:');
disp(t);
disp('2-stage gauss solution:');
disp(gauss2_sol);
disp('3-stage gauss solution:');
disp(gauss3_sol);
disp('exact:');
disp(y_exact);
%--------------------------------------------------
% plot
figure; hold on;
plot(t, y_exact, 'm-o', 'displayname','exact');
plot(t, gauss2_sol, 'r', 'displayname','2-stage gauss');
plot(t, gauss3_sol, 'g', 'displayname','3-stage gauss');
xlabel('t');
ylabel('y(t)');
title('Gauss Methods');
legend('location','best');
grid on;
hold off;

%% local function
function y = implicitRK_scalarLinear(a, b, c, t, f1, f2, y0)
% a- matrix of stage coeff ajm, b = weights of findal update
% cj vector nodes where stages rae evaluated
L = length(t) - 1;
y = zeros(1, L+1);
y(1) = y0;
s = length(b); % number of stages
for n = 1:L
    h = t(n+1) - t(n);
    M = zeros(s, s);
    r = zeros(s, 1);

    for j = 1:s % loopfor each stage
        alpha_j = f2(t(n) + c(j)*h);  %f2 at stage j, calcs how y affects the slope at the given stage, how much the current time step affected at t_n + cjh
%**************foreach row j, M(j,: is computed based on alpha coeff, step size h,
        for m = 1:s % inner loop
            if j==m
                M(j,m) = 1; % diag element =1 initially,
            end
            M(j,m) = M(j,m) - h*a(j,m)*alpha_j; %
        end
%***********
        r(j) = f1(t(n)+c(j)*h) + alpha_j*y(n); %sets the jth tney of the rhs
    end

    k = M\r; % solve s x s system
    sumStages = 0;
    for j = 1:s
        sumStages = sumStages + b(j)*k(j);
    end
    y(n+1) = y(n) + h*sumStages;
end
end
