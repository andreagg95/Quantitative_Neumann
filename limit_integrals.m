% This script evaluate the four terms in equation 3.27 in Gentile - Henrot
% and shows that the following inequalities hold
%
% \mu_1(h) \int_0^1 u_1(x)^2 h(x) dx < \mu_2(h) \int_0^1 u_2(x)^2 h(x)  dx
% \int_0^1 u_1(x)^2 h(x) dx > \int_0^1 u_2(x)^2 h(x) dx,
%

clc;
clear all;
close all;


% calculation of j01, omega2 and x0
f0 = @(x) besselj(0,x);
j01 = fzero(f0,2.4);
omega2 = 2*j01 + pi;
x0 = j01/omega2;


% definition of h
h = @(x) (x./x0).*(0 <= x & x < x0) + ...
         (x0 <= x & x <= 1-x0) + ...
         ((1-x)./x0).*(x> 1-x0 & x<=1);

%        ____________     
%      / |           |\
%     /  |           | \
%    /   |           |  \
%   /____|___________|___\
%   
%       x0         1-x0     
     

     
     
% computation of \mu1(h) solving the trascendental equation
%
% \tan((1-2*x0)*y) = (2*J_0(x_0*y)*J1(x_0*y))/(J_1(x_0*y)^2 - J_0(x_0*y)^2)
%

psi = @(y) (2*(besselj(0,x0.*y)).*(besselj(1,x0.*y)))/(besselj(1,x0.*y).^2 ...
            - besselj(0,x0.*y).^2);
psi2 = @(y) tan((1-2*x0).*y);
phi = @(y) psi2(y) - psi(y);
omega1 = fzero(phi,4.3);



% definition of u1 (3.26 Gentile-Henrot)
B = (besselj(0,omega1*x0))/(sin(omega1*x0-(omega1/2)));

u1 = @(x) besselj(0,omega1.*x).*(0 <= x & x < x0) + ...
         (B*sin(omega1.*x-omega1/2)).*(x0 <= x & x <=1-x0)  + ...
          - besselj(0,omega1.*(1-x)).*(1-x0 < x & x <=1);
     
     
     
% definition of u2 (3.26 Gentile-Henrot)  
u2 = @(x) besselj(0,omega2.*x).*(0 <= x & x < x0) + ...
         -(besselj(1,j01).*sin(omega2.*x-j01)).*(x0 <= x & x <=1-x0)+ ...
         + besselj(0,omega2.*(1-x)).*(1-x0 < x & x <=1);
     

% Normalization of u1 and u2 in L^2   
g1 = @(x) u1(x).^2;
g2 = @(x) u2(x).^2;

A1 = 1/(sqrt(integral(g1,0,1)));
A2 = 1/(sqrt(integral(g2,0,1)));


u1 = @(x) A1.*u1(x);
u2 = @(x) A2.*u2(x);



% computation of LHS1, RHS1, LHS2 and RHS2 (left and right hand side) that
% are 
%
%     LHS1 = \int_0^1 u1(x)^2*h(x) dx
%     RHS1 = \int_0^1 u2(x)^2*h(x) dx
%     LHS2 = mu1* \int_0^1 u1(x)^2*h(x) dx
%     RHS2 = mu2* \int_0^1 u2(x)^2*h(x) dx

% \mu_i = omega_i^2
mu1 = omega1^1;
mu2 = omega2^2;

% definition of the integrands
f1 = @(x) ((u1(x)).^2).*h(x);
f2 = @(x) ((u2(x)).^2).*h(x);

LHS1 = integral(f1,0,1);
RHS1 = integral(f2,0,1);

LHS2 = mu1*LHS1;
RHS2 = mu2*RHS1;


LHS1 > RHS1
LHS2 < RHS2

