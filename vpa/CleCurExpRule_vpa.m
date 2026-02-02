% CleCurExpRule_vpa - Standard and Modified Clenshaw-Curtis Rule
%
% VPA precision version
% 
% Computes the integral:
%
%   int_0^b f(s) exp(z s) ds
%
% Syntax:
%
%   integral = CleCurExpRule(fval, z, 0, b)
%
% ⚠ **THIS IS A NEW VERSION (2025) OF THE SUBROUTINE.** ⚠
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
% `fval` can also be an `(m+1) x n` matrix, where each column represents 
% function values for a different function. In this case, the function 
% returns `n` integrals, one for each column.
%
% -------------------------------------------------------------------------
% Error Estimation (for even m):
%
%   [integral, errorEst] = CleCurExpRule(fval, z, 0, b)
%
% If `fval` is an `(m+1) × n` matrix and `m` is **even**, the function also 
% returns `errorEst`, an estimate of the quadrature error. This is computed 
% by comparing results from the full set of nodes with those obtained 
% using the coarse mesh (`fval(1:2:end, :)`).
%
% **Limitations:**
% - If `m` is odd, `errorEst` is set to `NaN`, as error estimation is not 
%   supported in this case.
%
% This version employs MATLAB's Variable Precision Arithmetic (VPA) to
% enhance numerical stability when computing highly oscillatory integrals.
% As result, it is considerably slower than double precision implementation
% and it must be used only for testing purposes. 
% 
%
% -------------------------------------------------------------------------
% For more information:
% - https://www.arxiv.org/abs/2503.08169
%
% Author: Victor Dominguez
% Contact: victor.dominguez@unavarra.es
% Date: 2nd Feb 2025
% -------------------------------------------------------------------------
%
% Copyright (C) 2026 Victor Dominguez
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


function [integral, ErrEst]= CleCurExpRule_vpa(fval,z,varargin)

legacyUsed = ~isempty(varargin) && ...
             ~ischar(varargin{1}) && ...
             ~isstring(varargin{1});


op='CleCu'; %Clenshaw-Curtis rule by default


if legacyUsed
    % Legacy syntax requires b
    if numel(varargin) == 1
        % CleCurExpRule_vpa(fval, z, b)
        varargin = {'EndPoint', varargin{1}};

    elseif numel(varargin) == 2 && isnumeric(varargin{2}) && isscalar(varargin{2}) && varargin{2} == 2
        % CleCurExpRule_vpa(fval, z, b, 2)  -> Fejér rule
        varargin = {'EndPoint', varargin{1}, 'FejerRule', true};

    else
        error('CleCurExpRule_vpa:InvalidLegacySyntax', ...
              'Legacy syntax must be CleCurExpRule_vpa(f,z,b) or CleCurExpRule_vpa(f,z,b,2).');
    end
end

% ---------------------------------------------------------------------
% INPUT PARSER (VPA-safe)
% ---------------------------------------------------------------------
p = inputParser;

addParameter(p, 'EndPoint', sym(2), ...
    @(x) (isnumeric(x) || isa(x,'sym')) && isscalar(x));

addParameter(p, 'NumberOfNodes', 16, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 2 && mod(x,1) == 0);

addParameter(p, 'FejerRule', false, ...
    @(x) islogical(x) || isnumeric(x));

try
    parse(p, varargin{:});
catch ME
    error('CleCurExpRule_vpa:InvalidArgument', ...
          'Invalid input: %s', ME.message);
end

if legacyUsed
   % warning('CleCurExpRule_vpa:LegacySyntax', ...
   %     'Legacy positional syntax detected. Please use name-value arguments.');
end

% ---------------------------------------------------------------------
% Parsed parameters
% ---------------------------------------------------------------------
b      = vpa(p.Results.EndPoint);
m     = p.Results.NumberOfNodes;

 




if p.Results.FejerRule
    op='Fejer';
end 
if isa(fval, 'function_handle')
    if strcmp(op,'CleCu')
        xi = b/2 * (1 + cos((0:m) * pi/m));
        xi = xi(:);
    else
        xi = b/2*(1+cos((0.5:(m+0.5))*pi/(m+1)));
        xi = xi(:);
    end

    % Vectorized evaluation (preferred)
    try
        fval = fval(xi);
    catch
        % Non-vectorized 
        fval = arrayfun(fval, xi);
    end

elseif isvector(fval)
    fval = fval(:); % column

elseif ismatrix(fval)
    % OK

else
    error('CleCurExpRule:InvalidFval', ...
          'fval must be a vector, a (m+1)xN matrix, or a function handle.');
end



[m,n] = size(fval);
m = m-1; % m+1 nodes
znew = z*b/2;

if real(znew)>20
    warning('z in the reference interval [0,2] is greater than 20')
    warning('The integral is expected to exponentially large')
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
    if strcmp(op,'CleCu')
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
    fval = fval .* exp(znew*(xi+1));
   

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
integral = w.'*fval*b/2;
ErrEst = nan*integral;
if mod(m/2,1)==0 % We can compute an error estimate

    integral2 =  w2.'*fval(1:2:end,:)*b/2;
    ErrEst    = abs(integral-integral2);
end

return