% Segundo ensayo
% res = ensayo(z,alphas,n)
% 
% Compute 
% 
% \int_0^2 (s*(2-s))^alpha \exp(z s)\,{\rm d}s = 
%  double_factorial(2*alpha)*exp(z)*pi*besseli(alpha+1/2,z)*z^(-1/2-alpha)
%
% with the quad rule of n points
% Here alpha is 1/2, 3/2,... 
%
% Re z <= 0 (or positive and small)
%
% alpha can be a vector 
%
% Devuelve un ^
% 
% Feb 2025

function res = ensayo2(z,alphas,ns)
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
  
    t = linspace(0,pi,n+1);  t =1+ cos(t);
    y = g(t);
    [intN, ErrEst]= CleCurExpRule(y.',z,2);
    res{end}(end+1,1:4) = [n   intN  (intN-intex) abs(intN-intex)]  ;
    
   % disp(res(end,:))
    end
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
    
