%  [w, rho] = computeWeights(n, z)
%
% Computes the integrals:
%
%   w(n)    = int(T_{n-1}(s - 1), exp(z s),s,0,2)
%   rho(n)  = int(U_{n-1}(s - 1), exp(z s),s,0,2)
%
% The parameter `z` is assumed to satisfy Re(z) < r_0, where r_0 is small.
%
% [w, rho, condM] = CleCurRule(n, z)
%
% Returns `condM`, which represents the conditioning of the matrix 
% in the second phase of the algorithm.
%
% [w, rho, condM] = CleCurRule(n, z, n0)
%
% Uses the first-phase algorithm up to `n0`, and the second-phase 
% algorithm from `n0` onward.
%
% -------------------------------------------------------------------------
% For more information:
% - https://www.arxiv.org/abs/2503.08169
%
% Author:  Victor Dominguez
% Contact: victor.dominguez@unavarra.es
% Date:    11 March 2025
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



function [w,rho,condM] = computeWeights(n,z,n0)

condM = -1;
if nargin<3
    if real(z)==0
        n0 = max(ceil((abs(z)))+1,4);
    else
        n0 = 2*max(ceil(sqrt(abs(z)))+1,4);
        %disp(n0) 
    end
end

rho = zeros(n+1,1);

rho(1) = (exp(2*z)-1)/z;
rho(2) = (2*(1 +exp(2*z)*(-1 + z) + z))/z^2;

gamma      = 1/z*exp(2*z);
gamma2     = 1/z;

% First phase
for j = 2:min(n0,n+1)
    rho(j+1) = rho(j-1)-2*(j)/z*rho(j)+2*gamma +(-1)^(j+1)*2*gamma2;
end
% Second phase
if n>n0

    [rhoMax,nMax] = computeRhoMax(z,n);
    ns = sqrt(n0+1:nMax).';
    nsinv = 1./ns;
    nsinv =nsinv(2:end).*nsinv(1:end-1);

    Msize = length(ns)-1;

    tridiagCompressed=[ [-z/2*nsinv(2:end); 0] ...
        ones(length(nsinv),1) [0; z/2*nsinv(2:end)]];
    Matrix2 = spdiags(  tridiagCompressed,-1:1,Msize,Msize);

    rhs = zeros(Msize,1);
    rhs(1:2:end) = 2*gamma + (-1)^(n0+1)*2*gamma2;
    rhs(2:2:end) = 2*gamma - (-1)^(n0+1)*2*gamma2;
    rhs(1) = rhs(1)+rho(n0+1);
    rhs(end) = rhs(end)-rhoMax;
    rhs = ns(2:end).\rhs;


    if nargout ==3
        aux = Matrix2\rhs;
        condM = condest(Matrix2);
    else
        % Choose one:
        %aux = thomas_algorithm(tridiagCompressed(:,1), ...
        %            tridiagCompressed(:,2),tridiagCompressed(:,3), rhs);
        aux = Matrix2\rhs;
    end

    rho(n0+2:n0+Msize+1) =2/z*ns(2:end).\(aux);
end

% w weights
w    = zeros(n+1,1);
w(1) = rho(1);
w(2) = rho(2)/2;

w(3:2:n+1) =  -(2:2:n).'./z.*rho(2:2:n) +(gamma-gamma2);
w(4:2:n+1) =  -(3:2:n).'./z.*rho(3:2:n) +(gamma +gamma2);


rho = rho(1:n+1);
return



