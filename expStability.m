%% Built-in Function exp and expm1(x)

% Using exp(x), we can get at most ~ 55 digits for extremely small x. After
% that the result is not a good approximation of the real result. exp(x) is
% not stable because by the defination of stability, exp(x) fails to
% produce a result within machine epsilon when x <= 10^(-17), as shown
% below. And the error is not O(machine epsilon) as x becomes smaller and
% smaller. 

% expm1(x) is stable because it produces at least 100 digits. If we
% increase the accuracy, expm1(x) will give more accurate result than 100
% digits.

% In the Tylor expansion program for exp, despite that I use 1000 terms of
% expansion, the accurate result still vanishes when x <= 10^(-16), which
% has the same accuracy as the built-in exp(x). It is prbably due to the
% accuracy of other built-it functions that is called in EXP(x,n), such as
% '/' and factorial. 

%% Check the stability of exp(x)
close all; clc; clear;
n = 1;
y = 1; 
format long e
fprintf('  x                     exp(x)  \n')
% Let x range from 10^{-1} down to 10^{-20} %
while n <= 20   
    y = y/10;
    e = exp(y);
    n = n+1;
    fprintf('%13.16e %13.55e \n', y , e)
end

abs(exp(10^(-17)) - expm1(10^(-17))) / expm1(10^(-17)) 

%% Check the stability of expm1(x)
close all; clc; clear;
m = 1;
y = 1; 
format long e
fprintf('  x                     expm1(x)  \n')
% Let x range from 10^{-1} down to 10^{-30} %
while m <= 30
    y = y/10;
    e = expm1(y);
    m = m+1;
    fprintf('%13.16e %13.100e \n', y , e)
end

%% Use Tylor Expansion to Compute e^x
close all; clc; clear;

% Compute the Tylor expansion up to 1000 terms and produce the result for
% exp(x) for extremely small x. 

n = 1;
y = 1; 
format long e
fprintf('  x                     EXP(x,1000)  \n')
while n <= 20   
    y = y/10;
    e = EXP(y,1000);
    n = n+1;
    fprintf('%13.16e %13.55e \n', y , e)
end

% Input: 
% --- x the value we want to evaluate on
% --- n is the number of terms in the Taylor expansion 
function y = EXP(x, n)

format long

i = 0;
sum = 0;
while i <= n
    
    sum = sum + x^i/factorial(i);
    i = i + 1;
end 
 y = sum;
end













