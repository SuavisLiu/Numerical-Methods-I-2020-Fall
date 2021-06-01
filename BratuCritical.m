%% Problem 6
clear; clc; close all;
format short e
kmax = 200;
n = 128;
tol = 10^(-6);
lam = 3.51;

u1 = @(x) x.*(1-x);
t1 = linspace(0,1,n)';
x1 = highernewtons(u1, t1, n, lam, kmax, tol);
plot(t1,x1,'b')
hold on

u2 = @(x) 20*x.*(1-x);
t2 = linspace(0,1,n)';
x2 = highernewtons (u2, t2, n, lam, kmax, tol);
plot(t2,x2,'r')

xlabel('t')
ylabel('u(t)')
legend('u=t(1-t)', 'u=20*t(1-t)')
title({'Solution to 1d Bratu Problem with', ['lam = ', num2str(lam)]})


%% Higher Dimension Newton's Method
function x = highernewtons (u,t,n,lam,kmax,tol)
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
    v_p = [v(2:n);0];
    v_m = [0; v(1:n-1)];
    F = (1/h^2)*(-v_m + 2*v - v_p) - lam*exp(v);
    % compute the Jacobian at each iteration 
    D = 2 - lam*exp(v)*h^2;
    J = diag(D) + diag(-1*ones(1,n-1),1) + diag(-1*ones(1,n-1),-1);
    J = J/h^2;
    J = sparse(J);
    % x_k+1 by Newton's method
    v = v - J\F;
    fnorm(k+1) = norm(F);
    error = norm(v_prev - v);
    diffXk(k+1) = error;
    k = k+1;
end 
x = v;
NumOfIteration = [1:k]';
fnorm = fnorm';
diffXk = diffXk';
convRateOffNorm = [fnorm(2:k);0]./fnorm.^2.;
convRateOfX = [diffXk(2:k);0]./diffXk.^2.;
table(NumOfIteration,fnorm,convRateOffNorm)
table(NumOfIteration,diffXk,convRateOfX)
end
