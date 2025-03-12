% CleCurExpRule - Standard and Modified Clenshaw-Curtis Rule for 
%                 integral with exp(z s) term (z is complex)
%
% Computes the integral:
%
%   int(f(s) exp(z s),s,0,2)
%
% Syntax:
%
%   integral = CleCurExpRule(fval, z)
%
%     where `fval` is a COLUMN vector containing the function values at:
%
%         (1 + cos((0:m) * π / m))
%
%     These are the Chebyshev nodes.
%
%     ⚠ **The points MUST be provided as above, in REVERSE order.** ⚠
%
%   integral = CleCurExpRule(fval, z, 'EndPoint', b)
%
%     Computes the integral 
%
%         int(f(s) exp(z s),s,0,b)
%
%     where `b` is the endpoint, and `fval` is assumed to be a COLUMN vector 
%     containing the function values at:
%
%         (b/2) * (1 + cos((0:m) * π / m))
% 
%   integral = CleCurExpRule(f, z, 'EndPoint', b)
% 
%     Evaluates `f` at
%
%         (b/2) * (1 + cos((0:m) * π / m))
% 
%     and computes the integral over [0, b].
%
%   integral = CleCurExpRule(f, z, 'EndPoint', b, 'NumberOfNodes', m)
% 
%     Evaluates function `f` at `m+1` points given by:
%
%         (b/2) * (1 + cos((0:m) * π / m))
% 
%     and computes the integral in [0, b].  
%     ⚠ **Notice that the actual number of nodes is `NumberOfNodes = m + 1`.** ⚠
% 
% -------------------------------------------------------------------------
% Alternative Mode - Fejér Rule:
%
%   integral = CleCurExpRule(fval, z, 'FejerRule', 1,'EndPoint',b)
% 
% Computes the Fejér rule that evaluates the integral using the value
% of the function at 
%
%   (b/2) * (1 + cos((0.5:(m+0.5)) * π / (m+1)))
%
% and stored at fval
%
% ⚠ **FEJÉR RULE HAS NOT BEEN PROPERLY TESTED!** ⚠
%
% The Fejér rule does not include the endpoints of the interval as
% quadrature nodes, making it particularly suitable for functions with
% singularities at the interval’s extremes.
%
% -------------------------------------------------------------------------
% Matrix Input Support:
%
% `fval` can also be an `(m+1) × n` matrix, where each column represents
% function values for a different function. In this case, the function
% returns `n` integrals, one for each column.
%
% `f` can also be a ROW vector of `n` functions, in which case `n` integrals 
% are returned.
%
% -------------------------------------------------------------------------
% Error Estimation (for even m & Clenshaw-Curtis rule):
%
%   [integral, errorEst] = CleCurExpRule(fval, z)
%
% If `fval` is an `(m+1) × n` matrix and `m` is **even**, the function also
% returns `errorEst`, an estimate of the quadrature error. This is computed
% by comparing results from the full set of nodes with those obtained
% using the coarse mesh (`fval(1:2:end, :)`).
%
% **Limitations:**
% - If `m` is odd, `errorEst` is set to `NaN`, as error estimation is 
%   not yet supported in this case.
% -------------------------------------------------------------------------
% For more information:
% - https://www.arxiv.org/abs/2503.08169
%
% Author: Victor Dominguez
% Contact: victor.dominguez@unavarra.es
% Date: 11 March 2025
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

function [integral, ErrEst]= CleCurExpRule(fval,z,varargin)
op='CleCu'; %Clenshaw-Curtis rule

p = inputParser;
addParameter(p, 'EndPoint', 2, @(x) isreal(x) && x > 0);
addParameter(p, 'NumberOfNodes', 16, @(x)  isreal(x) && x > 1);
addParameter(p, 'FejerRule', false, @(x) isnumeric(x) || islogical(x)); 


try 
    parse(p, varargin{:});
catch ME
    errorParts = strsplit(ME.message, '.');
    errorMsg = errorParts{1}; % Get everything before the first period
    error('Custom Error: %s.\n', errorMsg);
end

b = p.Results.EndPoint;
if logical(p.Results.FejerRule)
    op='Fejer'; % Fejer Rule instead
end

if isa(fval,'function_handle')
   
    m = p.Results.NumberOfNodes;
    x= b/2 * (1 + cos((0:m) * pi/ m));
    x  = x(:); 
    fval = cell2mat(arrayfun(@(t) fval(t), x, 'UniformOutput', false));

elseif isvector(fval)
    fval = fval(:);
elseif ~ismatrix(fval)
    error('fval is neither a (m+1) x n matrix or a (vector of) funcion(s)')
end 

[m,n] = size(fval);
m = m-1;
znew = z*b/2;

if real(znew)>10
    warning('z in the reference interval [0,2] is greater than 10')
end

mEst=ceil(m/2);

if abs(znew)<1    % Classical rule

    if mod(m,2) == 0
        mEnd = m+1;
    else
        mEnd = m;
    end
    w = zeros(m+1,1);
    w(1) = 2;
    w(3:2:mEnd) = 2./(1-(2:2:mEnd).^2).';
    w2 = w(1:mEst+1);
    if op == 'CleCu'
        xi = cos((0:m)*pi/m).';
        w  = idctI(w);
        w2 = idctI(w2);
        % Correction for the first & last term
        w([1 end]) = 0.5*w([1 end]);
        w2([1 end]) = 0.5*w2([1 end]);
    else
        xi = cos((0.5:(m+0.5))*pi/(m+1)).';
        w  = idctII(w);
        w2 = idctII(w2);
    end
    fval = fval.* (exp(znew*(xi+1))*ones(1,n));
else
    w = computeWeights(m,znew);
    w2 = w(1:mEst+1);
    if op == 'CleCu'
        w  = idctI(w);
        w2 = idctI(w2);
        % Correction for the first & last term
        w([1 end])  = 0.5*w([1 end]);
        w2([1 end]) = 0.5*w2([1 end]);
    else
        w  = idctII(w);
        w2 = idctII(w2);
    end
end
integral = w.'*fval*b/2;
ErrEst = nan*integral;
if mod(m/2,1) ==0 & op == 'CleCu' % We can compute an error estimate
    integral2 = w2.'*fval(1:2:end,:)*b/2;
    ErrEst = abs(integral-integral2);
end

end