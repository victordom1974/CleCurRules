% CleCurExpRule - Standard and Modified Clenshaw-Curtis Rule
%
% VPA precision version
% 
% Computes the integral:
%
%   int_0^b f(s) exp(k s) ds
%
% Syntax:
%
%   integral = CleCurExpRule(fval, k, 0, b)
%
% ⚠ **THIS IS A NEW VERSION OF THE SUBROUTINE.** ⚠
%
% where `fval` is a COLUMN vector containing function values at:
%
%   b/2 * (1 + cos((0:m) * π / m))
%
% These are the Chebyshev nodes.  
%
% ⚠ **The points MUST be provided in REVERSE order.** ⚠
%
% -------------------------------------------------------------------------
% Alternative Mode - Fejér Rule:
%
%   integral = CleCurExpRule(fval, k, 0, 2)
%
% If the fourth argument is `2`, the **Fejér rule** is applied instead.
% In this case, `fval` must contain function values at:
%
%   b/2 * (1 + cos((0.5:(m+0.5)) * π / (m+1)))
%
% ⚠ **FEJÉR RULE IS CURRENTLY UNSUPPORTED!** ⚠
%
% The Fejér rule does not include the endpoints of the interval as 
% quadrature nodes, making it particularly suitable for functions with 
% singularities at the interval’s extremes.
%
% -------------------------------------------------------------------------
% Matrix Input Support:
%
% `fval` can also be an `n × (m+1)` matrix, where each column represents 
% function values for a different function. In this case, the function 
% returns `n` integrals, one for each column.
%
% -------------------------------------------------------------------------
% Error Estimation (for Odd m):
%
%   [integral, errorEst] = CleCurExpRule(fval, k, 0, b)
%
% If `fval` is an `m × n` matrix and `m` is **odd**, the function also 
% returns `errorEst`, an estimate of the quadrature error. This is computed 
% by comparing results from the full set of nodes with those obtained 
% using the coarse mesh (`fval(1:2:end, :)`).
%
% **Limitations:**
% - If `m` is even, `errorEst` is set to `NaN`, as error estimation is not 
%   supported in this case.
%
% This version employs MATLAB's Variable Precision Arithmetic (VPA) to
% enhance numerical stability when computing highly oscillatory integrals.
%
% -------------------------------------------------------------------------
% For more information:
% - https://www.arxiv.org/abs/2503.08169
%
% Author: Victor Dominguez
% Contact: victor.dominguez@unavarra.es
% Date: 05 March 2025
% -------------------------------------------------------------------------
%
% Copyright (C) 2025 Victor Dominguez
%
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% -------------------------------------------------------------------------


function [integral, ErrEst]= CleCurExpRule_vpa(fval,z,b,varargin)
op='CleCu'; %Clenshaw-Curtis rule


if nargin>3 && varargin{1}==2
    op='Fejer';
end
[m,n] = size(fval);
m = m-1;
integral = sym(zeros(1,n));
znew = z*b/2;

if real(znew)>20
    warning('k in the reference interval [0,2] is greater than 20')
end


mEst=ceil(m/2);


if abs(double(znew))<1    % Classical rule

    if mod(m,2)==0
        mEnd=m+1;
    else
        mEnd=m;
    end
   
    w = sym(zeros(m+1,1));
    w(1:2:mEnd)=2./(1-(0:2:mEnd).^2).';
    w = sym(w(:)); 
    w2 = w(1:mEst+1);
    if op == 'CleCu'
        xi = cos(sym((0:m)/m)*sym(pi)).';
        w  = idctI_vpa(w);
        w2 = idctI_vpa(w2);
        % Correction for the first & last term
        w([1 end]) = w([1 end])/sym(2);
        w2([1 end])= w2([1 end])/sym(2);
    else
        xi = cos(sym((0.5:(m+0.5))/(m+1))*sym(pi)).';
        w  = idctII_vpa(w);
        w2 = idctII_vpa(w2);
    end
    fval = fval.* (exp(znew*(xi+1))*ones(1,n));
   

else
    w = computeWeights_vpa(m,znew);
    %toc
    w2 = w(1:mEst+1);
    if op=='CleCu'
        w = idctI_vpa(w);
        w2= idctI_vpa(w2);
        % Correction for the first & last term
        w([1 end]) = w([1 end])/2;
        w2([1 end]) =w2([1 end])/2;
    else
        w  = idctII_vpa(w);
        w2 = idctII_vpa(w2);
    end


end 
integral= w.'*fval*b/2;
ErrEst=nan*integral;
if mod(m/2,1)==0 % We can compute an error estimate

    integral2= w2.'*fval(1:2:end,:)*b/2;
    ErrEst=abs(integral-integral2);
end

return