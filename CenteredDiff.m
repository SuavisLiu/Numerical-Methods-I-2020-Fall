%% Central Difference Quotient


% First, we modify h from 10^(-1) to 10^(-20) and compute the central
% difference quotient defined in Chapter 11. We obtain the error of the
% derivative by computing the absolute difference of the exact deriv and
% central difference quotient. The ratio in the table shows that when h is
% reduced by a factor of 10, the discretization error is reduced by a
% factor of about 100. The ratio is obtained by round the ratio of the
% previous error and the error under current h. 

% We see the ratio is appro. 100 when h is larger than 10^(-4). After that
% the error itself is not quite dominant by derivative values. The
% cancellation error begins to affect the result. When h becomes extremely
% small, i.e. h < 10^(-17), the central difference quotient is invisible
% due to the cancellation error and precision. 
% Therefore, approximately h ~ 10^(-4) gives the best results so that we
% are able to see both the central difference error and O(h^2) of the error. 


close all; clc; clear;
n = 1;
x = 1.0;
h = 1.0;
diffquo = 0.0;
error = 0.0;
deriv = cos(x);

fprintf(" deriv = %13.6e \n", deriv);
fprintf(" h      centraldiffquo  error         ratio\n");


% Let h range from 10^{-1} down to 10^{-20} %
while n <= 20
h = h/10; % h = 10^(-n) 
diffquo = (sin(x+h)-sin(x-h))/(2*h); % central difference quotient 
error = abs(deriv - diffquo);
if n ~= 1
    ratio = round(preverror/error);
else ratio = 0;
end
preverror = error;

fprintf("%5.1e %13.6e %13.6e   %d \n", h, diffquo, error, ratio);
n = n+1; 
end 
