%% Problem 6
% Helped by Bill Xu
%% part(a)
clear; close all; clc;
f = @(x) exp(-2*x)./(1+4*x);
a = 0; b = 1;
acc = integral(f,0,1); 
legpts = @(n) legpts(n,[a,b]); 
guassleg(f,legpts,acc,10)

%% part(b)
clear; close all; clc;
f = @(x) sin(x.^(1/3));
a = 10^(-3); b = 1;
acc = integral(f,a,b); 
legpts = @(n) legpts(n,[a,b]);
guassleg(f,legpts,acc,10);

%% part(c)
clear; close all; clc;
f = @(x) sin(x.^(1/3));
a = 10^(-6); b = 1;
acc = integral(f,a,b); 
legpts = @(n) legpts(n,[a,b]);
guassleg(f,legpts,acc,10);

%% part(d)
clear; close all; clc;
f = @(x) sin(x.^(1/3));
a = 0; b = 1;
acc = integral(f,a,b);
legpts = @(n) legpts(n,[a,b]);
guassleg(f,legpts,acc,10);


%% Guass Legendre Quadrature

function m = guassleg(f,legpts,int,tol) 
%
% Input: f -- the function; legpts -- the points xi's
%        int -- integral from Matlab Integral, to compare; tol -- 10^(-tol)
% Output: m -- the number of points we need to get 10^(-tol) accuracy
n_max = 1000;
m = 4;
for n = 1:n_max
    [nodes,weights] = legpts(m+1);
    integral = sum(f(nodes).*weights');
    if abs(integral-int) < 10^(-1*tol)
    break; 
    end
    m = m+1; 
end

end

