

% Fast discrete inverse cosine transform (IDCT-II)
% (Up to a multiplicative constant)
%
% Given (Y_k)_{k=0}^N Compute
%
% y_n=1/N*(0.5*Y(1)+\sum_{k=0}^{N-1} Y_k*cos(pi*(n+0.5)*k/N))
%
% Y is assumed to be a column vector with at least two
% elements
%
% Y can be a matrix. In this case the cosine transform
% is applied columnwise
%
% see http://en.wikipedia.org/wiki/Discrete_cosine_transform#DCT-II
%
% This function makes use of the native Matlab function idct 
%
% Last Update: 07 May 2013
function y=idctII_test(Y)
m = size(Y,1);
y2 = idct_fft_replacement(Y);
Y(2:end,:) = Y(2:end,:)*sqrt(2);
try 
    y = idct(Y);
    y = y/sqrt(m); 
    [y-y2]
catch
    y = y2;
end

end

function y = idct_fft_replacement(Y)
    % Reemplazo de IDCT con FFT para IDCT-II
    m = size(Y, 1);
    m0 = m-1;
   for j = 0:m0
        y(j+1, :) = 1/2 * Y(1, :) + ...
                     dot(Y(2:end, :)', cos(pi /(2* m) * (1:m0)' * (2*j+ 1)))';
   end
   y =y*2/m;

   W = (1:2:(2*m0+1))'*(1:(m0));
   W = cos(W*pi/(2*(m0+1)));
   y2 = 0.5*Y(1)+W*Y(2:end,:);
   y2 = y2*2/m;

   W2 = (0:(2*m0+1))'*(0:(2*m0+1))/(2*(m0+1));
   W2 = exp(1i*W2*pi);
   W2 = 1/m*W2;

  
   y3 = W2*[Y; 0*Y(1,:); -Y(end:-1:2,:)]; 
   y3 = y3(2:2:end,:);
   
   W4 = (0:(4*m0+3))'*(0:(4*m0+3))/(4*m0+4)*2;
   W4 = exp(1i*W4*pi)/m;
   Y4 = zeros(4*m,size(Y,2));
   Y4(1:m,:) = Y;

   Y4(end:-1:end-m+2,:) = Y(2:end,:);
   y4 = W4*Y4; 
   y4 = y4(2:2:2*m,:);

   y5 = fft(Y4)/m;
   y5 = y5(2:2:2*m,:); 
   y6 = mifft_vpa(Y4)/m;

   y6 = y6(2:2:2*m,:); 
   [y y2 y3 y4 y5 y6];
   y = y5;
   if isreal(Y)
       y = real(y);
   end
   
end