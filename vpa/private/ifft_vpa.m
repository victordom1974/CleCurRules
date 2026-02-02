% IFFT_VPA  Fast Inverse Fast Fourier Transform (IFFT)
%
% VPA precision version
%
%   Y = ifft_vpa(x)
%
% Computes the inverse discrete Fourier transform (IDFT) of a sequence
% using arbitrary precision arithmetic (VPA) or standard double precision,
% depending on the input type.
%
% Given a sequence $(x_n)_{n=0}^{N-1}$, this function computes:
%
%   Y_k = ∑_{n=0}^{N-1} x_n * exp(-2πi * n * k / N)
%
% (up to normalization constants depending on convention).
%
% The implementation supports both numeric and symbolic inputs.
%
%
% This implementation is mainly for educational and research purposes,
% not optimized for performance on large vectors.
%
% ### Reference
% * https://en.wikipedia.org/wiki/Fast_Fourier_transform
% * https://en.wikipedia.org/wiki/Discrete_Fourier_transform
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
% Date: 12 November 2025
function Y = ifft_vpa(x)


Y = x; 
   
N = size(x,1);

if mod(N,2)==1
    n = 0:N-1;    % Row vector for the time indices
    k = n;       % Column vector for the frequency indices

    % Create the DFT matrix using outer product
    if isa(x,'sym')

        W = vpa(exp(sym(-1i * 2 * pi / N * (n' * k))));
    else
        W = exp(-1i * 2 * pi / N * (n' * k));
    end
    % Compute the DFT as a matrix-vector multiplication
    Y = W * x;

else

    % Divide step
    even = ifft_vpa(x(1:2:end,:));  % Recur even indexed elements
    odd  = ifft_vpa(x(2:2:end,:));  % Recur odd indexed elements

    % Divide & Conquer step with the combination of even and odd part

    for k = 1:N/2
        % Complex exponential factor (twiddle factor)
        if isa(x,'sym')
            t = vpa(exp(sym(-2i*pi*(k-1)/N)));
        else
            t =  exp(-2i * pi * (k - 1) / N);
        end 
       
        Y(k,:) = even(k,:) + t* odd(k,:);
        Y(k + N/2,:) = even(k,:) - t*odd(k,:);
    end
 
end
