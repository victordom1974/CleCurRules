% First test
% 
% res = ensayo(z, ns)
%
% Returns in res:
%
% integral(@(t) pn(ns, t) .* exp(z * (t + 1)), -1, 1) 
%
% the approximation given by the quadrature rule with n+1 points and the error.
%
% IMPORTANT: Im(z) < 0
%
% The exact value is:
%
% (-1i)^(n) * (2 * pi / (1i * z))^(1/2) * exp(z) * besselj(n + 1/2, 1i * z)
%
% ns can be a vector of n values.
%
% Disclaimer: mostly in Spanish... 
%
% Assumptions:
% - Required files are available either in this directory or in the parent folder.
%
% For more information, refer to the Numerical Experiment section in:
% - https://www.arxiv.org/abs/2503.08169
%
% Feb 2025

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
%
% Mathematica command: 
%
% g2[n_, z_] = (-I)^n*(2*Pi/(I*z))^(1/2)*Exp[z]*BesselJ[n + 1/2, I*z]
res = []; 
for n = ns
    g = @(t) pn(n,t);

    intex = (-1i)^(n)*(2*pi/(1i*z))^(1/2)*exp(1i*imag(z))*besselj(n+1/2,1i*z,1);

    t = linspace(0,pi,n+1);  t = cos(t);
    y = g(t); 
    [intN, ErrEst]= CleCurExpRule(y.',z,2);
    res = [res; [intex intN abs(intex-intN)]]; 
   % disp(res(end,:))
   
end

return

function y = pn(n,x)
% Polinomio de legendre

y = legendre(n,x);
y = y(1,:);