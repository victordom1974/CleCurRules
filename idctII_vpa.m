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
# CleCurExpRule - Clenshaw-Curtis Quadrature with Exponential Weight

This MATLAB project implements Clenshaw-Curtis quadrature rules for computing integrals of the form:

$$ I = \int_0^b f(s) e^{z s} ds $$

## Installation
Clone the repository and add the files to your MATLAB path:

```matlab
addpath(genpath('CleCurExpRule'))
Functions Overview
Core Functions
These are the main functions implementing the quadrature rules:

CleCurExpRule – Computes the integral using the standard Clenshaw-Curtis rule.
computeWeights – Computes the quadrature weights for the Clenshaw-Curtis method.
computeRhoMax – Computes stability-related parameters.
Variable-Precision (VPA) Versions
These functions perform the same operations but using MATLAB’s Variable Precision Arithmetic (VPA) for higher accuracy:

CleCurExpRule_vpa – VPA version of CleCurExpRule.
computeWeights_vpa – VPA version of computeWeights.
computeRhoMax_vpa – VPA version of computeRhoMax.
Auxiliary Functions
These helper functions are used internally for efficient computation:

idctI – Implements the Type-I inverse discrete cosine transform.
idctII – Implements the Type-II inverse discrete cosine transform.
idctI_vpa – VPA version of idctI.
idctII_vpa – VPA version of idctII.
mifft_vpa – Performs the inverse FFT using variable-precision arithmetic.
Usage
To compute an integral using the standard Clenshaw-Curtis rule:

matlab
Copiar
Editar
result = CleCurExpRule(fval, z, b);
For higher precision using VPA:

matlab
Copiar
Editar
result_vpa = CleCurExpRule_vpa(fval, vpa(z), vpa(b));
License
This project is licensed under the MIT License - see the LICENSE file for details.

Author
Victor Dominguez
Email: victor.dominguez@unavarra.es
end

function y = idct_fft_replacement(Y)
% Reemplazo de IDCT con FFT para IDCT-II

m = size(Y, 1);

Y4 = sym(zeros(4*m,size(Y,2)));
Y4(1:m,:) = Y;
Y4(end:-1:end-m+2,:) = Y(2:end,:);

y = mifft_vpa(Y4)/m;

y = y(2:2:2*m,:);

if isreal(Y)
    y = real(y);
end

end