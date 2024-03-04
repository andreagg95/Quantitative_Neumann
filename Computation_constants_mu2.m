% this script evaluates the right hand side in 3.37 of Gentile, Henrot
% .....


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
     

     
     
% computation of \mu1(h) solving the trascendetal equation
%
% \tan((1-2*x0)*y) = (2*J_0(x_0*y)*J1(x_0*y))/(J_1(x_0*y)^2 - J_0(x_0*y)^2)
%

psi = @(y) (2*(besselj(0,x0.*y)).*(besselj(1,x0.*y)))/(besselj(1,x0.*y).^2 ...
            - besselj(0,x0.*y).^2);
psi2 = @(y) tan((1-2*x0).*y);
phi = @(y) psi2(y) - psi(y);
omega1 = fzero(phi,4.3);



% definition of mu1
B = (besselj(0,omega1*x0))/(sin(omega1*x0-(omega1/2)));

u1 = @(x) besselj(0,omega1.*x).*(0 <= x & x < x0) + ...
         (B*sin(omega1.*x-omega1/2)).*(x0 <= x & x <=1-x0)  + ...
          - besselj(0,omega1.*(1-x)).*(1-x0 < x & x <=1);
     
     
     
% definition of u2     
u2 = @(x) besselj(0,omega2.*x).*(0 <= x & x < x0) + ...
         -(besselj(1,j01).*sin(omega2.*x-j01)).*(x0 <= x & x <=1-x0)+ ...
         + besselj(0,omega2.*(1-x)).*(1-x0 < x & x <=1);
     

% Normalization of u1 and u2     
g1 = @(x) u1(x).^2;
g2 = @(x) u2(x).^2;

A1 = 1/(sqrt(integral(g1,0,1)));
A2 = 1/(sqrt(integral(g2,0,1)));

     
           

% definition of the integrands for C1 and C3
integrand1 = @(x) besselj(0,omega2.*x).*besselj(1,omega2.*x).*(x.^2);
integrand2 = @(x) besselj(0,omega2.*x).*besselj(1,omega1.*x).*(x.^2);


% definition of sigma
sigma = omega1/2 - omega1*x0;

% definition of the constants
C1 = (A2^2)*(omega2/x0)*integral(integrand1,0,x0);
C2 = 3*(A2^2)*((besselj(1,j01))^2)*(1-(2*x0))/2;
C3 = (A1*A2)*(omega1/x0)*integral(integrand2,0,x0);
C4 = (A1*A2)*(omega1*omega2)*(besselj(0,omega1*x0))*(besselj(1,j01))/((tan(sigma))*((omega2^2)-(omega1^2)));


% we write
% f(K2,c) = - (a20*K2^2 + a02*c^2 + a10*K2 + a01*c + a11*K2*c + a00)
% g(K2,c) = b20*K2^2 + b10*K2 + b11*K2*c

a20 = (12*C1 + 6*(C1/(x0^2)) -12*(C1/x0)+C2);
a02 = 12*(C1/(x0^2));
a10 = 6*(C1/(x0^2));
a01 = 12*(C1/(x0^2));
a11 = 12*(C1/(x0^2));
a00 = (4*C1)/(x0^2);

b20 = (6*(C3/(x0^2)) - 12*(C3/x0) + 6*C4);
b10 = (6*(C3/(x0^2)) - 12*(C3/x0) + 6*C4);
b11 = (12*(C3/(x0^2)) - 24*(C3/x0) + 12*C4);


% the function auxil is the characteristic function of the set
% { (K2,c) \in \R^2 : -1 <= c <= 0, -c-1 <= K2 <= -c }

aux = @(x,y) auxil(x,y);


% definition of the functions f and g
f1 = @(x,y) -( a20*(x.^2) + a02*(y.^2) + a10.*x + a01.*y + a11*(x.*y) + a00 );
g2 = @(x,y)  b20.*(x.^2) + b10.*x + b11*(x.*y);
g1 = @(x,y) abs(g2(x,y));

f = @(x,y) f1(x,y).*aux(x,y);
g = @(x,y) g1(x,y).*aux(x,y);


% definition of the difference
def = @(x,y) f(x,y) + g(x,y);


% f(K2,c) + |g(K2,c)| <= (max f(K2,c) + max |g(K2,c)|) = f(0,-1/2) + |g(1/2,-1)|

f(0,-1/2) + g(1/2,-1)


