%% -------------------------------------------------------------
%  Test script for CleCurExpRule
%  Different calling modes
%  + comparison with exact integral
%  + vector function implementation 
% --------------------------------------------------------------
%
% November 2025
%
% Last update: January 2026

clearvars   

%% PARAMETERS
f      = @(x) cos(6*x)-1i*sin(3*x);
z      = -14 + 25i;
b      = pi;
nNodes =32;  % Set nNodes even

%% Exact integral (reference)
syms x
ex_sym = int(f(x)*exp(z*x), x, 0, b);
ex     = double(ex_sym);

fprintf('Exact integral (reference): %.16e + %.16ei\n', real(ex), imag(ex));


%% --------------------------------------------------------------
% Method 1: Provide vector of function values manually
%% --------------------------------------------------------------

% Chebyshev nodes for interval [0, b]
t = (cos(pi*(0:nNodes)/nNodes) + 1).' * b/2;  % must be a column vector

ints = zeros(1,4);
errs = nan+zeros(1,4);

[ints(1), errs(1)] = CleCurExpRule(f(t), z, 'EndPoint', b);


%% --------------------------------------------------------------
% Method 2: Provide function handle, automatic node generation
%% --------------------------------------------------------------
[ints(2), errs(2)] =CleCurExpRule(f, z, ...
    'NumberOfNodes', nNodes, ...
    'EndPoint', b);


%% --------------------------------------------------------------
% Method 3: With error estimation
%% --------------------------------------------------------------
[ints(3), errs(3)] = CleCurExpRule(f, z, ...
    'NumberOfNodes', nNodes, ...
    'EndPoint', b);

%% --------------------------------------------------------------
% Method 4: With error estimation not available (because nNodes is odd)
%% --------------------------------------------------------------
[ints(4), errs(4)] = CleCurExpRule(f, z, ...
    'NumberOfNodes', nNodes+1, ...
    'EndPoint', b);





%% --------------------------------------------------------------
% RESULTS
%% --------------------------------------------------------------
fprintf('\n=== SCALAR TEST RESULTS ===\n');

signo = ['+','-'];
for k = 1:4
    fprintf('Method %d: %.16e %c %.16ei   (error est: %.3e)\n', ...
        k, real(ints(k)),  signo(1+double(imag(ints(k))>=0)),abs(imag(ints(k))), errs(k));
end

fprintf('\nExact:   %.16e + %.16ei\n\n', real(ex), imag(ex));


%% --------------------------------------------------------------
% Test with VECTOR-VALUED FUNCTION
%% --------------------------------------------------------------

clear errs_vec ints_vec
fvec = @(x) [5./(3+cos(3*x)),  exp(-2*cos(2*x))];

t = (cos(pi*(0:nNodes)/nNodes) + 1).' * b/2;

[ints_vec, errs_vec] = CleCurExpRule(fvec(t), z, ...
    'NumberOfNodes', nNodes, ...
    'EndPoint', b);


fvec = @(x) [5./(3+cos(3*x)),  exp(-2*cos(2*x))]; 

[ints_vec(end+1,:), errs_vec(end+1,:)] = CleCurExpRule(fvec(t), z, ...
    'NumberOfNodes', nNodes, ...
    'EndPoint', b);

% Fejer rules: no error estimate
[ints_vec(end+1,:), errs_vec(end+1,:)] = CleCurExpRule(fvec(t), z, ...
    'NumberOfNodes', nNodes, ...
    'EndPoint', b,'FejerRule',1);


[ints_vec(end+1,:), errs_vec(end+1,:)] = CleCurExpRule(fvec , z, ...
    'NumberOfNodes', nNodes, ...
    'EndPoint', b,'FejerRule',1);




fprintf('=== VECTOR FUNCTION TEST ===\n');
disp('Integrals:');
disp(ints_vec);
disp('Error estimates:');
disp(errs_vec);

