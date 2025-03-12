% THOMAS_ALGORITHM Solves a tridiagonal system Ax = d using the Thomas method.
%
% Inputs:
%   a - subdiagonal (vector of length n-1, shifted one element down)
%   b - main diagonal (vector of length n)
%   c - superdiagonal (vector of length n-1, shifted one element up)
%   d - vector of independent terms (of length n)
%
% Output:
%   x - solution of the tridiagonal system
%
% Example of usage:
%   a = [-1, -1, -1];       % subdiagonal (n-1 elements)
%   b = [2, 2, 2, 2];       % main diagonal (n elements)
%   c = [-1, -1, -1];       % superdiagonal (n-1 elements)
%   d = [1, 2, 3, 4];       % independent terms (n elements)
%   x = thomas_algorithm(a, b, c, d)
%
% This function supports symbolic variables, allowing computations
% in variable-precision arithmetic.
%
% GPT-4 writes most of this simple code, saving me a valuable
% amount of time.
% -------------------------------------------------------------------------
% For more information:
% - https://www.arxiv.org/abs/2503.08169
%
% Author: Victor Dominguez
% Contact: victor.dominguez@unavarra.es
% Date: 11 March 2025
% -------------------------------------------------------------------------
%
% Copyright (C) 2025 Victor Dominguez
%

function x = thomas_algorithm(a, b, c, d)

c_mod = c;
b_mod = b;
d_mod = d;

n = length(d); 
%matrix = spdiags([a(:) b(:) c(:)],-1:1,4,4);
%[l,u] =lu(matrix,0)
 
for i = 2:n
   % disp([i n])
    m = a(i-1) / b_mod(i-1);
    b_mod(i) = b_mod(i) - m * c(i);
    d_mod(i) = d_mod(i) - m * d_mod(i-1);
end

x = b;
x(n) = d_mod(n) / b_mod(n);

for i = n-1:-1:1

   % disp([i n])
    x(i) = (d_mod(i) - c_mod(i+1) * x(i+1)) / b_mod(i);
end

end
