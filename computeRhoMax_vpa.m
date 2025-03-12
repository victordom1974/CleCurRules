%  [rhoM, m] = computeRhoMax_vpa(z, varargin)
%
% VPA precision version
% 
% Computes the integral using Variable Precision Arithmetic (VPA):
%
%   rhoM = int_0^2 U_{m-1}(x - 1) * exp(z * x) dx
%
% where `m` is chosen to satisfy m > 1.6 * |z| for numerical stability.
%
% [rhoM, m] = computeRhoMax_vpa(z, m0)
%
% Ensures that `m` is greater than `m0`. The choice of `m` is guided by
% numerical stability considerations.
%
% [rhoM, m] = computeRhoMax_vpa(z, m0, M)
%
% Uses an `M x M` linear system to compute the approximation. If `M` is
% not specified, it is determined adaptively according to:
%
%   M = 2 * ceil(log(5 / (1e-32* abs(z) * r)) / log(1 + r)) + 2;
%
% where `r` is a parameter controlling the precision of the approximation
% (here, r = 2 is used).
%
% This version utilizes MATLAB's Variable Precision Arithmetic (VPA) to
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


function [rhoM,mc] = computeRhoMax_vpa(z,varargin)


 z=vpa(z,64); 
m = []; n0 = []; 

if nargin>1
    n0 = varargin{1};
    if nargin>2
        m  = varargin{2};
    end
end

r = 2;
if isempty(m)
    m  = 2*ceil(log(5/(1e-32*abs(z)*r))/log(1+r))+2; 
end

m1 = double(floor(max([1.5*abs(z)+2,2*8,n0+4])));
m2 = m1 + m; 
mc = (m1+m2)/2;
ms = sqrt(vpa(m1:m2)).';
msinv = 1./ms;
msinv = msinv(2:end).*msinv(1:end-1);

rhs = sym(zeros(double(m2-m1),1)); 

gamma      = 1/z*exp(2*z);
gamma2     = 1/z;
rhs(1:2:end) = 2*gamma + 2*gamma2;
rhs(2:2:end) = 2*gamma - 2*gamma2;

Msize = length(ms)-1;
tridiagCompressed=[ [-z/2*msinv(2:end); 0] ...
            ones(length(msinv),1) [0; z/2*msinv(2:end)]]; 


%Matrix2 = spdiags([ [-z/2*msinv(2:end); 0] ones(length(msinv),1) z/2*msinv],-1:1,Msize,Msize);





rhs = (ms(2:end).\rhs); 
%y0 = (Matrix2\rhs); 
y0 =  thomas_algorithm(tridiagCompressed(:,1), tridiagCompressed(:,2),tridiagCompressed(:,3), rhs); 

rhoM = 1/sqrt(vpa(mc))*z/2*y0((m2-m1)/2);
return
%e1 = zeros(m2-m1,1);
%e2=e1;
%e1(1)=1;
%e2(end)=1;
%r=Matrix2\[e1 e2];
%

% rhs(1:2:end) = 2*gamma + (-1)^(n0+1)*2*gamma2;
%    rhs(2:2:end) = 2*gamma - (-1)^(n0+1)*2*gamma2;
%
%
%r((m2-m1)/2,:)
%n = m1+(m2-m1)/2