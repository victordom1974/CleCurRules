
z = -250;
ns = 10:1:500;

computeQ = 1;

if computeQ

    zsExperiment = z*exp(-1i*pi/2*[1 4/5 3/5 2/5 1/5 0]);

    % Ensure that zs is is truly real or truly imaginary
    % (eliminating round-off effects).
    % ind1 = abs(real(zsExperiment)./abs(zsExperiment))<1e-14;
    % zsExperiment(ind1) =  -1i*abs(imag(zsExperiment(ind1)));

    % Fix to enforce a negative imaginary part when z is purely imaginary.
    % This is required for computing the "true integral" consistently,
    % due to the branch of the Bessel function in the exact integral.
    % ind2 = abs(imag(zsExperiment)./abs(zsExperiment))<1e-14;
    % zsExperiment(ind2) = real(zsExperiment(ind2));

    results = [];
    for zs = zsExperiment

        disp("z = " + num2str(zs))
        res = ensayo(zs,ns);
        res = round(real(res),3,'significant')+...
            1i*round(imag(res),3,'significant');

        results{end+1} = [ns',res];
    end
end

figure(1)
clf
subplot(321)
figure(1)
titles{1} = "$z = " + num2str(z) +"{\rm i}$"
titles{2} = "$z = " + num2str(z) +"\exp( {\rm i}\pi/5) $"
titles{3} = "$z = " + num2str(z) +"\exp( 2{\rm i}\pi/5)$ "
titles{4} = "$z = " + num2str(z) +"\exp( 3{\rm i}\pi/5) $"
titles{5} = "$z = " + num2str(z) +"\exp( 4{\rm i}\pi/5)$ "
titles{6} = "$z = " + num2str(z) + "$ "
sgtitle('Exact & Computed integral (real & imaginary part)' );

for j =1:length(results)
    subplot(3,2,j)
    title(titles{j})
    res = results{j}
    ns = res(:,1);
    p = plot(ns, real(res(:,2)),ns,imag(res(:,2)),...
        ns, real(res(:,3)),ns,imag(res(:,3)),'LineWidth',1)
    %set(gca,'yscale','log')
    axis tight

    title(titles{j},'Interpreter','latex','FontSize',14)
end

figure(2)
sgtitle('Modulus of the (absolute) error' );
for j =1:length(results)
    subplot(3,2,j)
    title(titles{j})
    res = results{j}
    ns = res(:,1);
    p = plot(ns, abs(res(:,4)),'LineWidth',1)
    %set(gca,'yscale','log')
    axis tight

    title(titles{j},'Interpreter','latex','FontSize',14)
end


figure(3)

sgtitle('Relative error' );
for j =1:length(results)
    subplot(3,2,j)
    title(titles{j})
    res = results{j};
    ns = res(:,1);
    res = round(abs(res(:,4)./res(:,2)),3,'significant');
    p = plot(ns, res,'LineWidth',1)
    set(gca,'yscale','log')
    axis tight
    ylim([1e-18,1])
    xlim([0,ns(end)/1.5])
    title(titles{j},'Interpreter','latex','FontSize',14)

end