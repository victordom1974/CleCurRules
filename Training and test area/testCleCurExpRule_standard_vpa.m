% Check the Product Clenshaw-Curtis rule package 
%
% This script evaluates the integral:
%
%   ∫_0^2 f(s) exp(z s) ds
%
% where z can be a complex number.
%
% It is implicitly assumed that Re(z) is small or moderate, 
% since even moderate values of Re(z) can make the integral extremely large.
%
% Two implementations of the rule are provided:
%
% - Double precision.
% - VPA precision (32-digit precision) for testing purposes.
%
% The second implementation requires the MATLAB Symbolic Toolbox. 
%
% -------------------------------------------------------------------------
% Copyright (C) 2025 Victor Dominguez
%
% This script is provided "as is" without any express or implied warranties
%
% You may use, modify, and distribute this script for academic and 
% research purposes with proper attribution.
%
% Author: Victor Dominguez
% Contact: victor.dominguez@unavarra.es
% Date: March 2025
% -------------------------------------------------------------------------


% Insert your favourite function. 
% 
% 

% Store the current directory
originalDir = pwd;
% Change to the parent directory and add it to the path
parentDir = fullfile(originalDir, '..');
addpath(parentDir); 

symbolic = 0
% Parameters: z 
z = -500+95i;  
val =[];

f = @(x) cos(4*x).*exp(-x);
if symbolic == 0
    syms x
    ex = int(f(x)*exp(z*x),x,0,2);
else
    ex = []; % We will use the result of the rule in the finest mesh
             % as exact value
end

% number of nodes for the experiment:
ns = [8 16 32 64 128 256 512 1024];
ns_vpa = [ 8 16 32]; % For omiting ns_vpa=[]

% Chose one rule:
% Clenshaw
%Cle = 1; disp('Clenshaw-Curtis rule')
% Fejer
Cle = 0; disp('Fejer rule')


% Exact integral.
% Need the symbolic toolbox 

for n = ns
    if Cle == 1
        t= linspace(0,pi,n+1).';
        t = cos(t)+1;

        y = f(t);
        res = CleCurExpRule(y,z);
    else
        t= linspace(0,pi,n+2).';
        t = (t(1:end-1)+t(2:end))/2;
        t = cos(t)+1;

        y = f(t);
        res = CleCurExpRule(y,z,'FejerRule',1);
    end
    val = [val; res];
end
if ~isnan(ex)
    error = val-ex;
    nExperiments = size(val,1);
else
    error = val(1:end-1,:)-val(end,:);
    nExperiments = size(val,1)-1;
end
disp('Results')
formattedStr = '';
for k = 1:nExperiments 
    if imag(val(k)) >= 0
        formattedStr = sprintf('%s N = %5.3d  %24.16e + %24.16ei\n', formattedStr,ns(k), real(val(k)), abs(imag(val(k))));
    else
        formattedStr = sprintf('%s N = %5.3d  %24.16e - %24.16ei\n', formattedStr,ns(k), real(val(k)), abs(imag(val(k))));
    end
end
disp(formattedStr);

disp('Error')
formattedStr = '';
for k = 1:nExperiments 
    if imag(val(k)) >= 0
        formattedStr = sprintf('%s N = %5.3d  %24.16e + %24.16ei\n', formattedStr,ns(k), real(error(k)), imag(error(k)));
    else
        formattedStr = sprintf('%s N = %5.3d %24.16e - %24.16ei\n', formattedStr, ns(k), real(error(k)),  imag(error(k)));
    end
end
disp(formattedStr);

disp('Version VPA')
disp('===========')
disp(' ')
val=[];
for n =ns_vpa
    if Cle == 1
        t = linspace(sym(0), sym(pi), sym(n) + 1).';
        t = cos(t)+1;
        y = vpa(f(t),64); 

        res = [CleCurExpRule_vpa(y,vpa(z),2)];
    else
        t = linspace(sym(0), sym(pi), sym(n) + 2).';
        t = (t(1:end-1)+t(2:end))/2;
        t = cos(t)+1;

        y = vpa(f(t),64);
        res = CleCurExpRule_vpa(y,vpa(z),2,2);
    end



    val = [val; res];
end
if ~isnan(ex)
    error = val-ex;
    nExperiments = size(val,1);
else
    error = val(1:end-1,:)-val(end,:);
    nExperiments = size(val,1)-1;
end
disp('Results')

formattedStr = '';
for k = 1:nExperiments 
    if imag(val(k)) >= 0
        formattedStr = sprintf('%s N = %5.3d  %24.16e + %24.16ei\n', formattedStr,ns_vpa(k), real(error(k)), imag(error(k)));
    else
        formattedStr = sprintf('%s N = %5.3d %24.16e - %24.16ei\n', formattedStr, ns_vpa(k), real(error(k)), abs(imag(error(k))));
    end
end 

disp(formattedStr);
disp('Error')
formattedStr = '';
for k = 1:nExperiments 
    if imag(val(k)) >= 0
        formattedStr = sprintf('%s N = %5.3d  %24.16e + %24.16ei\n', formattedStr,ns_vpa(k), real(error(k)), imag(error(k)));
    else
        formattedStr = sprintf('%s N = %5.3d %24.16e - %24.16ei\n', formattedStr, ns_vpa(k), real(error(k)), abs(imag(error(k))));
    end
end
disp(formattedStr);