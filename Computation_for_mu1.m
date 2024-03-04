% Function to evaluate the integral in 3.14 in Gentile-Henrot
%
% int_0^{1/2} ( J_1(\omega x)^2 - J_0(\omega x)^2 ) x^3 dx


clc;
clear all;
close all;


% define J_0 and j_{0,1}
f0 = @(x) besselj(0,x);
j01 = fzero(f0,2.4);

% definition of \omega = sqrt(mu1)
omega = 2*j01;


% definition of the integrand
g = @(x) (((besselj(1,omega.*x)).^2)- ((besselj(0,omega.*x)).^2)).*(x.^3);

% evaluation
(omega^2)*integral(g,0,1/2)




