%% Problem 5  Symmetry 
clear; close all; clc
format long e
kmax = 100;
n = 128;
tol = 10^(-6);

u1 = @(x) x.*(1-x);
t1 = linspace(0,1,n+2)';
t1 = t1(2:end-1); 
t1 = t1(1:n/2);
T = linspace(0,1,n)';
x1 = highernewtons (u1, t1, n, kmax, tol);
X1 = [x1;fliplr(x1')']; % mirror the second half
plot(T, X1,'b')
hold on 

u2 = @(x) 20*x.*(1-x);
x2 = highernewtons (u2, t1, n, kmax, tol);
X2 = [x2;fliplr(x2')']; % mirror the second half
plot(T,X2,'r')


xlabel('t')
ylabel('u(t)')
legend('u=t(1-t) with k = 3', 'u=20*t(1-t) with k = 6')
title('Solution to 1d Bratu Problem Exploit Symmetry')

%% Higher Dimension Newton's Method with Symmetry 
function x = highernewtons (u, t, n, kmax, tol)
%
% Inputs: u is the initial guess in form of function handle; 
%         t is the initial set of points
%         n is the number of intervals in [0,1]
%         kmax is the max of iteration 
%         tol is the tolerant error
% Output: x = u(t_k) is solution to -u" - exp(u) = 0
% 
h = 1/(n+1); % step size
% Initial guess v_0= u(t_0)
v = u(t);
k = 0; % number of iteration 
error = 10;
while k < kmax && error > tol
    v_prev = v;
    % compute F at each iteration    
    v_p = [v(2:end);v(end)]; %64
    v_m = [0; v(1:end-1)]; 
    F = (1/h^2)*(-v_p + 2*v - v_m) - exp(v);
    % compute the Jacobian at each iteration 
    D = 2 - exp(v)*h^2;
    D = [D(1:end-1); 1 - exp(v(end))*h^2];
    J = diag(D) + diag(-1*ones(1,n/2-1),1) + diag(-1*ones(1,n/2-1),-1);
    J = J/h^2;
    J = sparse(J);
    % x_k+1 by Newton's method
    v = v - J\F;
    error = norm(v_prev - v);
    k = k+1;
end 
x = v;
k
end
