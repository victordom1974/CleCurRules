% Fast Inverse Discrete Cosine Transform (iDCT-I)
%
% VPA precision version
%
% Given a sequence (y_k)_{k=0}^{N}, this function computes:
%
%   Y_k = (1/N) * (0.5 * (y_0 + (-1)^k * y_N) + ...
%          \sum_{n=1}^{N-1} y_n * cos(π * n * k / N))
%
% The input `y` is assumed to be a column vector with at least two elements.
%
% If `Y` is a matrix, the transform is applied column-wise.
%
% Reference: https://en.wikipedia.org/wiki/Discrete_cosine_transform#DCT-I
%
% -------------------------------------------------------------------------
% Copyright (C) 2025 Victor Dominguez
%
% This function is provided "as is" without any express or implied warranties.
%
% You may use, modify, and distribute this script for academic and 
% research purposes with proper attribution.
%
% Author: Victor Dominguez
% Contact: victor.dominguez@unavarra.es
% Date: March 2025
% -------------------------------------------------------------------------

function Y=idctI_vpa(y)

if size(y,1)==1
    Y=y;
    y=y(:);
else
    Y=y;
end
m=size(y,1);

y=[y; y(end-1:-1:2,:)];
z=mifft_vpa(y)/(m-1);

% Test if y is real. If so, consider only the real part of FFT
ind = ~any(imag(double(y)));
if ind
    z(:,ind)=real(z(:,ind)); 
end
Y(:)=z(1:m,:);
