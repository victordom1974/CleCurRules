
z0 = -40;
ns = [80 160 320 640 1280 2560 5120];
rs =0:4;
alphas = [1/2 3/2 5/2];
for j=0:3
    clear results

    results = cell(1, length(alphas));
    for k = 1:length(alphas)
        results{k} = ns';
    end
    for r = rs
        zr= z0*exp(-1i*j*pi/6)*4^(r);
        aux = ensayo2c(zr,alphas, ns);
        for l = 1:length(results)
            results{l}(:,end+1) = aux{l}(2:end,end);
        end
    end
    disp(' ')
    disp(" |z|="+abs(z0)+"*4^r"+ " exp(-1i*" + num2str(j)+ "*pi/6)")
    disp(' ')


    disp('e.o.c')
    disp('=====')
    disp(' ')
    for j= 1: length(alphas)
        disp("alpha = " + alphas(j))
        tab = results{j};
        printtab(tab);
        log2(results{j}(1:end-1,2:end)./results{j}(2:end,2:end))
    end
    %printtab(tab);
    disp("=============")

end
