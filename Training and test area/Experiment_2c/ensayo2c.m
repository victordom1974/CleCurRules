% ENSAYO2C
%
%   res = ensayo2c(z, alpha, ns)
%
% Compute the integral
%
%     ∫_0^2 (s(2 - s))^α exp(z s) ds
%
% whose exact value is
%
%     double_factorial(2α) * exp(z) * pi * besseli(α + 1/2, z) * z^(-1/2 - α).
%
% The approximation is obtained using an (ns+1)-point product Cle–Cur
% quadrature rule.
% Valid values of alpha are 1/2, 3/2, 5/2, ...
%
% INPUT:
%   z      : complex parameter in the exponential
%   alpha  : scalar or vector of admissible exponents (1/2, 3/2, 5/2, …)
%   ns     : number(s) of quadrature points (actually ns+1 nodes are used)
%
% OUTPUT:
%   res    : if alpha is a scalar, res is a matrix with the structure
%              res(1,:)   = [alpha, z, exact_value]
%              res(2:end,:) = [ns; computed_integrals; errors; abs(errors)]
%            if alpha is a vector, res is a cell array with length(alpha),
%            and res{j} contains the matrix corresponding to alpha(j).
%
% Last update: November 2025



function res = ensayo2c(z,alphas,ns)
p = cd;
cd ..
addpath(cd);
cd(p)


res = {};
for alpha = alphas
    g = @(t) (t.*(2-t)).^alpha ;

    intex =double_factorial(2*alpha)*pi*exp(imag(z)*1i)...
        *besseli(alpha+1/2,z,1)*z^(-1/2-alpha);
    res{end+1}= [alpha z intex];

    for n = ns

        t    = linspace(0,pi,n+1);  t =1+ cos(t);
        y    = g(t);
        intN = CleCurExpRule(y.',z);
        res{end}(end+1,1:4) = [n   intN  (intN-intex) abs(intN-intex)]  ;

        % disp(res(end,:))
    end
end

if length(res)==1
    res = res{1};
end
return

function y = pn(n,x)
% Polinomio de legendre

y = legendre(n,x);
y = y(1,:);


function df = double_factorial(n)
if n == 0 || n == 1
    df = 1;
else
    df = prod(n:-2:1); % Producto de n, (n-2), (n-4), ..., 1 o 2
end

