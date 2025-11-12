% -------------------------------------------------------------------------
% CleCurExpRule - Standard and Modified Clenshaw–Curtits Rule 
%                 for integrals with exp(z s) weight
%                 (double precision version)
% -------------------------------------------------------------------------
%
% Computes the integral
%
%       ∫_0^b f(s) exp(z s) ds
%
% in a robust and flexible way using:
%   - Standard Clenshaw–Curtis rule
%   - Modified exponential-weight Clenshaw–Curtis variant
%   - Optional Fejér rule (experimental, not fully validated)
%
% -------------------------------------------------------------------------
% BASIC USAGE
%
%   integral = CleCurExpRule(fval, z)
%
% where `fval` is a COLUMN vector containing function values at Chebyshev
% nodes in REVERSE order:
%
%       (1 + cos((0:m) * π / m))
%
% These nodes correspond to the interval [0,2].
%
% -------------------------------------------------------------------------
% INTEGRAL IN A GENERIC INTERVAL [0, b]
%
%   integral = CleCurExpRule(fval, z, 'EndPoint', b)
%
% where `fval` contains the values at
%
%       (b/2) * (1 + cos((0:m) * π / m))
%
% -------------------------------------------------------------------------
% AUTOMATIC NODE GENERATION
%
%   integral = CleCurExpRule(f, z, 'EndPoint', b, 'NumberOfNodes', m)
%
% Evaluates function handle `f` at the Chebyshev nodes and computes:
%
%       ∫_0^b f(s) exp(z s) ds
%
% ⚠ The actual number of quadrature points is m+1.
%
% -------------------------------------------------------------------------
% FEJÉR RULE (experimental)
%
%   integral = CleCurExpRule(fval, z, 'FejerRule', true, 'EndPoint', b)
%
% Uses nodes:
%
%       (b/2) * (1 + cos((0.5:(m+0.5)) * π / (m+1)))
%
% ⚠ FEJÉR RULE HAS NOT BEEN PROPERLY TESTED.
%
% -------------------------------------------------------------------------
% MATRIX INPUT SUPPORT
%
% - `fval` can be (m+1) × n → returns n integrals.
% - `f` can be a ROW vector of n function handles → returns n integrals.
%
% -------------------------------------------------------------------------
% ERROR ESTIMATION (Clenshaw–Curtis only, even m)
%
%   [integral, errorEst] = CleCurExpRule(fval, z)
%
% If:
%   - the rule is Clenshaw–Curtis (not Fejér)
%   - m is even
%
% then:
%
%       errorEst = |I_full – I_coarse|
%
% where the coarse mesh uses every second node.
%
% Otherwise, `errorEst = NaN`.
%
% -------------------------------------------------------------------------
% SEE ALSO
%   https://www.arxiv.org/abs/2503.08169
%
% Author:  Victor Dominguez
% Contact: victor.dominguez@unavarra.es
% Date:    11 March 2025
% -------------------------------------------------------------------------
%
% Copyright (C) 2025
% Permission is hereby granted under the MIT License.
% -------------------------------------------------------------------------


function [integral, ErrEst] = CleCurExpRule(fval, z, varargin)

% ---------------------------------------------------------------------
% INPUT PARSER
% ---------------------------------------------------------------------
p = inputParser;

addParameter(p, 'EndPoint', 2, ...
    @(x) isnumeric(x) && isscalar(x) && x > 0);

addParameter(p, 'NumberOfNodes', 16, ...
    @(x) isnumeric(x) && isscalar(x) && x >= 2 && mod(x,1)==0);

addParameter(p, 'FejerRule', false, ...
    @(x) islogical(x) || isnumeric(x));

try
    parse(p, varargin{:});
catch ME
    error('CleCurExpRule:InvalidArgument', ...
          'Invalid input: %s', ME.message);
end

b = p.Results.EndPoint;
m = p.Results.NumberOfNodes;

% Default rule
op = 'CleCu';
if logical(p.Results.FejerRule)
    op = 'Fejer';
end


% ---------------------------------------------------------------------
% HANDLE FUNCTION INPUT (fval is a function handle)
% ---------------------------------------------------------------------
if isa(fval, 'function_handle')

    x = b/2 * (1 + cos((0:m) * pi/m));
    x = x(:);

    % Vectorized evaluation (preferred)
    try
        fval = fval(x);
    catch
        % Non-vectorized 
        fval = arrayfun(fval, x);
    end

elseif isvector(fval)
    fval = fval(:); % column

elseif ismatrix(fval)
    % OK

else
    error('CleCurExpRule:InvalidFval', ...
          'fval must be a vector, a (m+1)xN matrix, or a function handle.');
end


% ---------------------------------------------------------------------
% DIMENSIONS AND z SCALING
% ---------------------------------------------------------------------
[m, n] = size(fval);
m = m - 1;                         % number of subintervals
if m < 1
    error('CleCurExpRule:TooFewNodes', ...
          'NumberOfNodes must be >= 2.');
end

znew = z * b/2;

if real(znew) > 20
    warning('CleCurExpRule:LargeZ', ...
            'z in reference interval [0,2] exceeds 20, integral may be large.');
end

mHalf = ceil(m/2);


% ---------------------------------------------------------------------
% SMALL |z| → CLASSICAL RULE
% ---------------------------------------------------------------------
if abs(znew) < 1

    % Define highest odd index
    if mod(m,2)==0
        mEnd = m+1;
    else
        mEnd = m;
    end

    w = zeros(m+1,1);
    w(1) = 2;
    w(3:2:mEnd) = 2./(1 - (2:2:mEnd).^2);

    w2 = w(1:mHalf+1);

    if strcmp(op,'CleCu')
        xi = cos((0:m)*pi/m).';

        w  = idctI(w);
        w2 = idctI(w2);

        % endpoint corrections
        w([1 end])  = 0.5 * w([1 end]);
        w2([1 end]) = 0.5 * w2([1 end]);

    else   % Fejér
        xi = cos((0.5:(m+0.5))*pi/(m+1)).';
        w  = idctII(w);
        w2 = idctII(w2);
    end

    fval = fval .* exp(znew*(xi+1));

% ---------------------------------------------------------------------
% LARGE |z| → MODIFIED WEIGHTS
% ---------------------------------------------------------------------
else
    w  = computeWeights(m, znew);
    w2 = w(1:mHalf+1);

    if strcmp(op,'CleCu')
        w  = idctI(w);
        w2 = idctI(w2);

        w([1 end])  = 0.5 * w([1 end]);
        w2([1 end]) = 0.5 * w2([1 end]);

    else
        w  = idctII(w);
        w2 = idctII(w2);
    end
end


% ---------------------------------------------------------------------
% FINAL INTEGRAL
% ---------------------------------------------------------------------
integral = (w.' * fval) * b/2;

% ---------------------------------------------------------------------
% ERROR ESTIMATION
% ---------------------------------------------------------------------
ErrEst = NaN * integral;

if strcmp(op,'CleCu') && mod(m,2)==0
    integral2 = (w2.' * fval(1:2:end,:)) * b/2;
    ErrEst    =  integral - integral2;
end

end
