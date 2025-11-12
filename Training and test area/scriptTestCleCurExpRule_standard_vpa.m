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


% Insert your "favourite" escalar or vectorial (row) function. 
% 
% 
f = @(x) [cos(4*x).*exp(-x)  (3+sin(5*x))];
% f = @(x) [cos(4*x)    (3+sin(5*x))];

% Store the current directory
originalDir = pwd;
% Change to the parent directory and add it to the path
parentDir = fullfile(originalDir, '..');
addpath(parentDir); 

%symbolic = 1; % The exact integral is computed analytically 
symbolic = 0; %The "exact" integral is computed numerically 
% Parameters: z 
z = -500+95i;  
val =[];

if symbolic == 1
    syms x
    ex = int(f(x)*exp(z*x),x,0,2);
else
    ex = []; % We will use the result of the rule in the finest mesh
             % as exact value
end

% number of nodes for the experiment:
ns = 1+[8 16 32 64 128 256 512 1024];
ns_vpa = [ 8 16 32 64 ]*15/8; % For omiting ns_vpa=[]
ns_vpa = [ 15 31 63 127]+1;
% Chose one rule:
% Clenshaw
Cle = 1; disp('Clenshaw-Curtis rule')
% Fejer
%Cle = 0; disp('Fejer rule')


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
    error = val-double(ex);
    nExperiments = size(val,1);
else
    error = val(1:end-1,:)-val(end,:);
    nExperiments = size(val,1)-1;
end
disp('Results')
formattedStr = '';
for k = 1:nExperiments
    formattedStr = sprintf('%s N = %5.3d  |',formattedStr,ns(k));
    for ell = (1:size(error,2)) 
        re = real(val(k,ell));
        im = imag(val(k,ell));
        signo = ['-','+'];
        formattedStr = ...
          sprintf('%s  %c%18.16e %c %18.16ei  |', formattedStr,...
          signo((re<0)+1), abs(re), signo((im<0)+1), abs(im));
    end
    formattedStr = sprintf('%s \n',formattedStr);

end  
disp(formattedStr);

disp('Error')
formattedStr = '';
for k = 1:nExperiments
    formattedStr = sprintf('%s N = %5.3d   |',formattedStr,ns(k));
    for ell = (1:size(error,2)) 
        re = real(error(k,ell));
        im = imag(error(k,ell));
        signo = ['-','+'];
        formattedStr = ...
          sprintf('%s  %c%4.3e %c %4.3ei  |', formattedStr,...
          signo((re<0)+1), abs(re), signo((im<0)+1), abs(im));
    end
    formattedStr = sprintf('%s \n',formattedStr);

end  
disp(formattedStr);

disp('Version VPA')
disp('===========')
disp(' ')
val=[];
for n =ns_vpa
    disp(n)
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
    formattedStr = sprintf('%s N = %5.3d  |',formattedStr,ns(k));
    for ell = (1:size(error,2)) 
        re = real(val(k,ell));
        im = imag(val(k,ell));
        signo = ['-','+'];
        formattedStr = ...
          sprintf('%s  %c%18.16e %c %18.16ei  |', formattedStr,...
          signo((re<0)+1), abs(re), signo((im<0)+1), abs(im));
    end
    formattedStr = sprintf('%s \n',formattedStr);

end  
disp(formattedStr);
 
disp('Error')
formattedStr = '';
for k = 1:nExperiments
    formattedStr = sprintf('%s N = %5.3d   |',formattedStr,ns(k));
    for ell = (1:size(error,2)) 
        re = real(error(k,ell));
        im = imag(error(k,ell));
        signo = ['-','+'];
        formattedStr = ...
          sprintf('%s  %c%4.3e %c %4.3ei  |', formattedStr,...
          signo((double(re)<0)+1), abs(re), signo((double(im)<0)+1), abs(im));
    end
    formattedStr = sprintf('%s \n',formattedStr);

end  
disp(formattedStr);