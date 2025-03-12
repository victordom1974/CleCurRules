% Stability Test for Computing Weights rho(z)
%
% This script tests the numerical stability of the computation of 
% the weights rho(z), from which w(z) is derived.
%
% Author: Victor Dominguez
% Date: March 2025
%
% Assumptions:
% - Required files are available either in this directory or 
%   in the parent folder.
%
% For more information, refer to the Numerical Experiment section in:
% - https://www.arxiv.org/abs/2503.08169
%

% We compute the following:
%
%   rho_stab:        Computed using the fully stable algorithm.
%   rho_vpa_stab:    Computed using the fully stable algorithm 
%                    with quadruple precision. Required the symbolic
%                    toolbox
%   rho_nonstab:     Computed using the first algorithm, which is unstable 
%                    for n > |z|^{1/2} or n > |z| if z is purely imaginary.
%   rho_vpa_nonstab: Computed using the first algorithm with quadruple 
%                    precision. It has the same instability as the previous
%                    one.
%
% Next, we generate graphs illustrating these behaviors.



% Store the current directory
originalDir = pwd;
% Change to the parent directory and add it to the path
parentDir = fullfile(originalDir, '..');
addpath(parentDir); 

rehacer =0
if rehacer
    n = 256;
    zs0 = -40;
    zs =   [0 1/12 1/6 1/4];
    clf
    rho_nonstab = [];
    rho_stab = [];
    rho_vpa_stab = [];
    rho_vpa_nonstab =[];

    w_nonstab = [];
    w_stab = [];
    w_vpa_stab = [];
    w_vpa_nonstab =[]

    for j = 1:length(zs)
        z = zs0*exp(2*pi*1i*sym(zs(j)));
        z = z *sym(pi);
       % z = double(z);
       n0 = double(2*max(ceil(sqrt(abs(z)))+1,4));
        [w,rho] = computeWeights(n,double(z),n0);
        w_stab{end+1} = w;
        rho_stab{end+1} = rho;
        [w,rho] = computeWeights(n,double(z),n);
        w_nonstab{end+1} = w;
        rho_nonstab{end+1} = rho;

        [w,rho] = computeWeights_vpa(n,double(z));
        w_vpa_stab{end+1} = w;
        rho_vpa_stab{end+1} = rho;
        [w,rho] = computeWeights_vpa(n,double(z),n);
        w_vpa_nonstab{end+1} = w;
        rho_vpa_nonstab{end+1} = rho;
    end
end

% Plotting
figure(1)
titles{1} = "$z = - 40\pi $"
titles{2} = "$z = - 40\pi\exp( {\rm i}\pi/3) $"
titles{3} = "$z = - 40\pi\exp( 2{\rm i}\pi/3)$ "
titles{4} = "$z = - 40\pi \rm i$ "
for j = 1:4

    subplot(2,2,j)
    % Bien
    %set(gca,'yscale','linear')
    set(gca,'yscale','log')
    axis tight
    ns = 0:length(rho_nonstab{j})-1;
    p = plot(ns,abs(rho_nonstab{j}),ns,abs(rho_vpa_nonstab{j}));

    ylim([1e-5,1e11])
    xlim([0,n])
    yticks= [ 1e-4 1e-1, 1e02, 1e5,1e8, 1e11] ;
    set(gca,'ytick',yticks)
    set(gca,'yscale','log')

    legend("Double prec.", 'Quadruple prec.','fontsize',10,'location','northwest')
    title(titles{j},'Interpreter','latex','FontSize',14)
    ax = gca;
    grid on
    drawnow
    ax.YMinorGrid = 'off';
end


% Plotting
figure(2)
titles{1} = "$z = - 40\pi $"
titles{2} = "$z = - 40\pi\exp( {\rm i}\pi/3) $"
titles{3} = "$z = - 40\pi\exp( 2{\rm i}\pi/3)$ "
titles{4} = "$z = - 40\pi {\rm i}$ "
for j = 1:4

    subplot(2,2,j)
    % mal
    %set(gca,'yscale','linear')
    set(gca,'yscale','log')
    axis tight

    p = plot(ns,abs(rho_vpa_stab{j}-rho_nonstab{j}),ns,abs(rho_vpa_stab{j}-rho_vpa_nonstab{j}));
    ydata = get(p,'Ydata');
    for r=1:length(ydata)
        ydata{r} = round(ydata{r},3,'significant');
    end
    %set(p,'ydata',ydata{end})

    ylim([1e-32,1e8])
    xlim([0,n])
    %    set(gca,'ytick',[1e-10, 1e-06, 1e-02, 1e02, 1e06])
    yticks= [1e-28, 1e-22,1e-16,1e-10,1e-4,1e2,1e8];
    set(gca,'yscale','log')
    set(gca,'ytick',yticks)


    legend("Double prec.", 'Quadruple prec.','fontsize',10,'location','northwest')
    title(titles{j},'Interpreter','latex','FontSize',14)
    ax = gca;
    grid on
    drawnow
    ax.YMinorGrid = 'off';
end


% Plotting
figure(3)
titles{1} = "$z = - 40\pi $"
titles{2} = "$z = - 40\pi\exp( {\rm i}\pi/3) $"
titles{3} = "$z = - 40\pi\exp( 2{\rm i}\pi/3)$ "
titles{4} = "$z = - 40\pi {\rm i}$ "
for j = 1:4

    subplot(2,2,j)
    % mal
    %set(gca,'yscale','linear')
    set(gca,'yscale','log')
    axis tight

    p = plot(ns,abs(rho_vpa_stab{j}-rho_stab{j}));
    ylim([1e-19,1e-14])
    xlim([0,n])
    %    set(gca,'ytick',[1e-10, 1e-06, 1e-02, 1e02, 1e06])
    yticks= [1e-18,1e-17,1e-16, 1e-15];
    set(gca,'yscale','log')
    set(gca,'ytick',yticks)


    title(titles{j},'Interpreter','latex','FontSize',14)
    ax = gca;
    grid on
    drawnow
    ax.YMinorGrid = 'off';
    legend('Error rho\_vpa\_stab-rho\_stab')
end

% Plotting
figure(4)
titles{1} = "$z = - 40\pi $"
titles{2} = "$z = - 40\pi\exp( {\rm i}\pi/3) $"
titles{3} = "$z = - 40\pi\exp( 2{\rm i}\pi/3)$ "
titles{4} = "$z = - 40\pi \rm i$ "
for j = 1:4

    subplot(2,2,j)
    % Bien
    %set(gca,'yscale','linear')
    set(gca,'yscale','log')
    axis tight
    ns = 0:length(rho_nonstab{j})-1;
    p = plot(ns,abs(rho_stab{j} -rho_vpa_nonstab{j}));

    ylim([1e-18,1e-010])
    xlim([0,ns(end)])
    yticks= [ 1e-16 1e-14 1e-12 1e-10] ;
    set(gca,'ytick',yticks)
    set(gca,'yscale','log')

    %legend("Error between Algorithm 0,'fontsize',10,'location','northwest')
    title(titles{j},'Interpreter','latex','FontSize',14)
    ax = gca;
    grid on
    drawnow
    ax.YMinorGrid = 'off';

    legend('Error rho\_vpa\_stab-rho\_vap\_nonstab')
end
