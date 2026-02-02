% IDCTII  Fast Inverse Discrete Cosine Transform (IDCT-II)
%
% (Up to a multiplicative constant)
%
% VPA precision version
%
% y = idctII_vpa(Y)
% 
% Given a sequence (Y_k)_{k=0}^{N}, this function computes:
%
%   y_n = (1/N) * (0.5 * Y_0 + ∑_{k=0}^{N-1} Y_k * cos(π * (n + 0.5) * k / N))
%
% The input `Y` is assumed to be a column vector with at least two elements.
%
% If `Y` is a matrix, the transform is applied column-wise.
%
% Reference: https://en.wikipedia.org/wiki/Discrete_cosine_transform#DCT-II
%
% This function utilizes MATLAB's built-in `idct` function.
%
% The user can choose between:
% - MATLAB's native implementation.
% - An alternative approach based on FFT, which can be extended to 
%   symbolic variables.
%
% -------------------------------------------------------------------------
% Copyright (C) 2025 Victor Dominguez
%
% This function is provided "as is" without any express or implied
% warranties.
%
% You may use, modify, and distribute this script for academic and 
% research purposes with proper attribution.
%
% Author: Victor Dominguez
% Contact: victor.dominguez@unavarra.es
% Date: 05 March 2025
% -------------------------------------------------------------------------


function y = idctII_vpa(Y)
m = size(Y,1);

y = idct_fft_replacement(Y);

end

function y = idct_fft_replacement(Y)
% Reemplazo de IDCT con FFT para IDCT-II

m = size(Y, 1);

Y4 = sym(zeros(4*m,size(Y,2)));
Y4(1:m,:) = Y;
Y4(end:-1:end-m+2,:) = Y(2:end,:);

y = ifft_vpa(Y4)/m;

y = y(2:2:2*m,:);

if isreal(Y)
    y = real(y);
end

end
