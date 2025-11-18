% Uncomment to add this to the path
% addpath(genpath(fileparts(fileparts(pwd))))

%Compute the integral for a smooth case 
f = @(x) cos(5*pi*x)./(4+sin(4*pi*x));

% zs are given below
%
% This corresponds to the second example in case
ns = [10 20 40 80 160 320 640 2560];

res = 0;
clear integ results
js = 0:3
for j = js

    zs = -20*[1 4 16 64 256 1024] * exp(pi*1i*j/6);
    r = 1; 
    
    for n = ns 
        s = 1; 
        for z = zs
            
            %t = linspace(0,1,n+1);
            %x = cos(pi*t);
            %fval = f(x+1);
            %integval = CleCurExpRule(fval(:), z);

            integval = CleCurExpRule(f, z, 'NumberOfNodes', n);
      
            % Fejér rule:
            % integval = CleCurExpRule(@(x) f(x+1), z, 'NumberOfNodes', n, 'FejerRule', 1);

            integ(r,s) = integval;
            s = s+1; 

        end
        r = r+1; 
    end

    results{j+1} = abs(integ(1:end-1,1:end) - integ(end,:));

    M = round(results{end}, 3, 'significant'); % Round matrix to significant digits

    disp("j= "+j)
    fprintf('\\begin{bmatrix} \n'); % Start LaTeX matrix environment
    fprintf('n \\ z &')
    fprintf('\n')
    for i = 1:size(M,1)
        fprintf('%3d   & %5.2e', ns(i), M(i,1)); % Print first element without "&"
        for j = 2:size(M,2)
            fprintf('  &   %5.2e', M(i,j)); % Separate values using "&"
        end
        fprintf(' \\\\\n'); % New LaTeX row
    end
    fprintf('\\end{bmatrix}\n'); % End LaTeX matrix environment
end
