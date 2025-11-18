% ENSAYO  Compute quadrature approximation and error for a model integral.
%
%   res = ensayo(z, ns) returns the quadrature approximation of
%
%       ∫_{-1}^{1} p_n(s) * exp( z * (s + 1) ) ds
%
%   together with the corresponding error. p_n is the Lagrange polynomial.
%  
%   The computation uses a quadrature rule with n+1 nodes, where n is 
%   each entry of NS.
%
%   INPUT:
%     z:    Complex parameter with Im(Z) < 0.
%     ns:   Scalar or vector of non-negative integers. Each value n in NS
%           specifies the degree of p_n and the number of quadrature points
%           (n+1) to be used.
%
%   OUTPUT:
%     res(k,1)  - Exact integral.
%     res(k,2)  - Computed quadrature approximation.
%     res(k,3)  - Module of the error .
%
%   The exact value of the integral is
%
%       (-1i)^n * ( 2*pi/(1i*z) )^(1/2) * exp(z) * besselj(n+1/2, 1i*z).

%
%   Example:
%       z  = -10 - 5i;
%       ns = 1:5;
%       res = ensayo(z, ns);
%
%   Nov 2025

function res = ensayo(z,ns)
p = cd;
cd ..
addpath(cd);
cd(p)

%res = integral(@(t)  pn(n,t).*exp(z*(t+1)),-1,1) - ...
%    (1i)^(n)*(2*pi/(1i*z))^(1/2)*exp(z)*besselj(n+1/2,-1i*z);

% for n =1:10
%  [ integral(@(t)  pn(n,t).*exp(z*(t+1)),-1,1) ,
%     (-1i)^(n)*(2*pi/(1i*z))^(1/2)*exp(z)*besselj(n+1/2,1i*z)]
% end
% for n =1:10
%  [ integral(@(t)  pn(n,t).*exp(z*(t+1)),-1,1) ,
%     (-1i)^(n)*(2*pi/(1i*z))^(1/2)*exp(1i*imag(z))*besselj(n+1/2,1i*z,1)]
% % end
% g2[n_, z_] = (-I)^n*(2*Pi/(I*z))^(1/2)*Exp[z]*BesselJ[n + 1/2, I*z]
res = []; 
for n = ns
    g = @(t) pn(n,t);

    intex = (-1i)^(n)*(2*pi/(1i*z))^(1/2)*exp(1i*imag(z))*besselj(n+1/2,1i*z,1);

    t = linspace(0,pi,n+1);  t = cos(t);
    y = g(t);
    intN = CleCurExpRule(y.',z,'endpoint',2);
    res = [res; [intex intN abs(intex-intN)]]; 
   % disp(res(end,:))
   
end

return

function y = pn(n,x)
% Polinomio de legendre

y = legendre(n,x);
y = y(1,:);