%% -------------------------------------------------------------
%  Test script for CleCurExpRule
%  Different calling modes 
%  + Fejer vs CleCur product rule
%  + comparison with exact integra
%  + Vector function implementation
% --------------------------------------------------------------
%
% January 2026
%

clearvars   

%% PARAMETERS
f      = @(x) cos(6*x)-1i*sin(3*x);
z      = 0.3;
b      = pi;
nNodes = 20;  % Set nNodes even for error estimates for the CleCur rule

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

% Chebyshev-fejer nodes for interval [0, b]
 
t2 = b/2*(1+cos((0.5:(nNodes+0.5))*pi/(nNodes+1))).';

intCleCur = zeros(1,4);
intFejer = intCleCur;
errs = nan+zeros(1,4);

intCleCur(1) = CleCurExpRule(f(t), z, 'EndPoint', b);


%% --------------------------------------------------------------
% Method 2: Provide function handle, automatic node generation
%% --------------------------------------------------------------
intCleCur(2) = CleCurExpRule(f, z, ...
    'NumberOfNodes', nNodes, ...
    'EndPoint', b);

% Method 2a: Fejer rule is selected
%% --------------------------------------------------------------
intFejer(1) =  CleCurExpRule(f(t2), z, ...
    'NumberOfNodes', nNodes, ...
    'EndPoint', b,'FejerRule',1);

intFejer(2) =  CleCurExpRule(f, z, ...
    'NumberOfNodes', nNodes, ...
    'EndPoint', b,'FejerRule',true);

%% --------------------------------------------------------------
% Method 3: With error estimation
%% --------------------------------------------------------------
[intCleCur(3), errs(3)] = CleCurExpRule(f, z, ...
    'NumberOfNodes', nNodes, ...
    'EndPoint', b);

%% --------------------------------------------------------------
% Method 4: With error estimation not available (because nNodes is odd)
%% --------------------------------------------------------------
[intCleCur(4), errs(4)] = CleCurExpRule(f, z, ...
    'NumberOfNodes', nNodes+1, ...
    'EndPoint', b);



%% --------------------------------------------------------------
% RESULTS
%% --------------------------------------------------------------
fprintf('\n=== SCALAR TEST RESULTS ===\n');

signo = ['+','-'];
fprintf('Clenshaw-Curtis product rule\n')
for k = 1:4
    fprintf('Method %d: %.16e %c %.16ei   (error est: %.3e)\n', ...
        k, real(intCleCur(k)),  signo(1+double(imag(intCleCur(k))>=0)),abs(imag(intCleCur(k))), errs(k));
end


fprintf('Fejer product rule\n')

for k = 1:2
    fprintf('Method %d: %.16e %c %.16ei  \n', ...
        k, real(intFejer(k)),  signo(1+double(imag(intFejer(k))>=0)),abs(imag(intFejer(k))));
end
fprintf('\nExact:   %.16e + %.16ei\n\n', real(ex), imag(ex));


%% --------------------------------------------------------------
% Test with VECTOR-VALUED FUNCTION
%% --------------------------------------------------------------

clear errs_vec intCleCur_vec intFejer_vec errs_vec


fvec = @(x) [5./(3+cos(3*x)),  exp(-2*cos(2*x))]; 


% Chebyshev nodes for interval [0, b]
t = (cos(pi*(0:nNodes)/nNodes) + 1).' * b/2;  % must be a column vector

% Chebyshev-fejer nodes for interval [0, b]
t2 = (0:(nNodes+1))/(nNodes+1);
t2 = pi/2*(t2(2:end)+t2(1:end-1));
t2 = (cos(t2)+ 1).' * b/2;  % must be a column vector


[intCleCur_vec, errs_vec] = CleCurExpRule(fvec(t), z, ...
    'NumberOfNodes', nNodes, ...
    'EndPoint', b);


intFejer_vec = CleCurExpRule(fvec(t2), z, ...
    'NumberOfNodes', nNodes, ...
    'EndPoint', b,'FejerRule',true);

intCleCur_vec(end+1,:) = CleCurExpRule(fvec, z, ...
    'NumberOfNodes', nNodes, ...
    'EndPoint', b);

intFejer_vec(end+1,:) = CleCurExpRule(fvec , z, ...
    'NumberOfNodes', nNodes, ...
    'EndPoint', b,'FejerRule',1);



disp(" ")


fprintf('=== VECTOR FUNCTION TEST ===\n');
disp('Integrals with Exponential Product Clenshaw-Curtis rule:');
disp(intCleCur_vec);

disp('Integrals with Exponential Product Fejer rule:');
disp(intFejer_vec);
disp('Error estimates:');
disp(' - Exponential Product Clenshaw-Curtis rule')
disp(errs_vec);

