
z0 = -40;
ns = [80 160 320 640 1280 2560 5120]; 
rs =0:4;
alphas = [3/2 3/2 ];
for j=0:3
    clear results
    results{1} = ns';
    results{2} = ns'; 
    for r = rs
        zr= z0*exp(-1i*j*pi/6)*4^(r);
        aux = ensayo2(zr,alphas, ns);
        aux{1}
        for l = 1:length(results)
            results{l}(:,end+1) = aux{l}(2:end,end);
        end
    end
    disp(' ')
    disp(" |z|="+abs(z0)+"*4^r"+ "z= 4^r*("+z0+")")
    disp("alpha = "+alphas(1))
    tab = results{1};
    printtab(tab);
    disp("alpha = "+alphas(2))
    tab= results{2};
    disp('e.o.c')
    log2(results{2}(1:end-1,:)./results{2}(2:end,:))
    printtab(tab);
    disp("=============")
    pause
    
end
