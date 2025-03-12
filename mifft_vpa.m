function Y = mifft_vpa(x)


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
    Y(:) = W * x;

else

    % Divide step
    even = mifft_vpa(x(1:2:end,:));  % Recur even indexed elements
    odd  = mifft_vpa(x(2:2:end,:));  % Recur odd indexed elements

    % Conquer step with the combination of even and odd
    combined = x*0;
    for k = 1:N/2
        % Complex exponential factor (twiddle factor)
        if isa(x,'sym')
            t = vpa(exp(-2i * pi * (k - 1) / N));
        else
            t =  exp(-2i * pi * (k - 1) / N);
        end
        Y(k,:) = even(k,:) + t* odd(k,:);
        Y(k + N/2,:) = even(k,:) - t*odd(k,:);
    end
 
end
