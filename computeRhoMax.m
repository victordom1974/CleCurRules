%  [rhoM, m] = computeRhoMax(z, m0)
%
% Computes the integral:
%
%   rhoM = int_0^2 U_{m-1}(x - 1) * exp(z * x) dx
%
% where `m` is chosen to satisfy n > 1.8 * |z|.
%
% [rhoM, m] = computeRhoMax(z, m0)
%
% Ensures that `m` is greater than `n0`. The choice of `m` is guided by
% numerical stability considerations.
%
% [rhoM, m] = computeRhoMax(z, m0, M)
%
% Uses an `M x M` linear system to compute the approximation. If `M` is
% not specified, it is determined adaptively according to:
%
%   M = 2 * ceil(log(5 / (1e-16 * abs(z) * r)) / log(1 + r)) + 2;
%
% where `r` is a parameter controlling the precision of the approximation.
% (r = 1.8 is taken here)
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

function [rhoM,mc] = computeRhoMax(z,varargin)

 
m = []; m0 = []; 

if nargin>1
    m0 = varargin{1};
    if nargin>2
        m  = varargin{2};
    end
end

r = 1.8;
if isempty(m)
    m  = 2*ceil(log(5/(1e-16*abs(z)*r))/log(1+r))+2; 
end


m1 = double(floor(max([1.2*abs(z)+2,2*8,m0+4])));
m2 = m1 +m; 
ms = sqrt(m1:m2).';
msinv = 1./ms;
msinv =msinv(2:end).*msinv(1:end-1);

rhs = zeros(m2-m1,1); 

gamma      = 1/z*exp(2*z);
gamma2     = 1/z;
rhs(1:2:end) = 2*gamma + 2*gamma2;
rhs(2:2:end) = 2*gamma - 2*gamma2;

Msize = length(ms)-1;

tridiagCompressed=[ [-z/2*msinv(2:end); 0] ...
            ones(length(msinv),1) [0; z/2*msinv(2:end)]]; 

rhs = (ms(2:end).\rhs);
y0 = thomas_algorithm(tridiagCompressed(:,1), tridiagCompressed(:,2),tridiagCompressed(:,3), rhs);
 

mc = (m1+m2)/2;
rhoM = z/2*1./sqrt(mc)*y0((m2-m1)/2);
% Alternative 
% Matrix2 = spdiags(  tridiagCompressed,-1:1,Msize,Msize);  
% y2 = z/2* (ms(2:end).\(Matrix2\rhs) ); 
% rhoM = y2((m2-m1)/2)

end