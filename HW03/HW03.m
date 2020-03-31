%% AE 370 HW3
% Linyi Hou (linyih2)

close all hidden
clear;clc

nvect = [5,10,25,50];
% nvect = 3;

fcn = @(x) 1./(1+x.^2);

err = zeros(length(nvect),1);

for j = 1 : length(nvect)
    
    n = nvect(j);
    
    xx = linspace(-5,5,1000)';
    grid_ones = ones(size(xx));
    
    b = zeros(n+1,1);
    
    % DEFINITIONS
    % phi  - orthogonal basis used to define G
    % b    - R.H.S. of the linear system Gc = b
    
    %Build first two rows of b and first two terms of polynomial
    %approximant fa
    %(slightly different for k = 1 & 2...)
    phi_1h = grid_ones ./ trapz(xx,grid_ones.*grid_ones);
    phi_1  = phi_1h ./ sqrt(trapz(xx,phi_1h.^2));
    b(1) = trapz(xx,fcn(xx).*phi_1);
    fa = b(1) .* phi_1;
    
    phi_2h = xx - trapz(xx,xx.*grid_ones)/trapz(xx,grid_ones.*grid_ones);
    phi_2  = phi_2h ./ sqrt(trapz(xx,phi_2h.^2));
    b(2) = trapz(xx,fcn(xx).*phi_2);
    fa = fa + b(2) .* phi_2;
    
    %remaining n-1 rows
    phi_km2 = phi_1;
    phi_km1 = phi_2;
    for jj = 3 : n+1
        
        phi_kh = xx .* phi_km1 - ...
               trapz(xx,xx.*phi_km1.*phi_km1)./trapz(xx,phi_km1.^2).*phi_km1 - ...
               trapz(xx,xx.*phi_km1.*phi_km2)./trapz(xx,phi_km2.^2).*phi_km2;
        
        phi_k = phi_kh ./ sqrt(trapz(xx,phi_kh.^2));
        
        b(jj) = trapz(xx,fcn(xx).*phi_k);
        
        fa = fa + b(jj) .* phi_k;
        
        phi_km2 = phi_km1;
        phi_km1 = phi_k;
    end
    
    %Compute error
    err(j) = norm( fcn(xx) - fa );
    
    figure(1)
    subplot(2,2,j)
    plot( xx, fcn(xx), '-', 'linewidth', 2 ), hold on
    plot( xx, fa, '--', 'linewidth', 2 )
    ylim([-0.2 1.2])

    %make plot pretty
    title( ['$n = ', num2str( n ),'$'] ,'interpreter', 'latex',...
        'fontsize', 16)
    if j == 1
        h = legend( '$f(x)$', '$f_a(x)$');
    end

    if j <= 2
        set( gca, 'XTick', [] )
    else
        xlabel( '$x$', 'interpreter', 'latex', 'fontsize', 16)

    end
    set(h, 'location', 'NorthEast', 'Interpreter', 'Latex', 'fontsize', 16 )
    set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )

    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'Color', [1 1 1])
    set(gca, 'Color', [1 1 1])
    set(gcf, 'PaperUnits', 'centimeters')
    set(gcf, 'PaperSize', [25 25])
    set(gcf, 'Units', 'centimeters' )
    set(gcf, 'Position', [0 0 25 25])
    set(gcf, 'PaperPosition', [0 0 25 25])
    
    grid(gca,'minor')
    grid on
    
end

figure(1)
print( '-dpng', 'vary_n', '-r200' )

%plot error
figure(100)
semilogy( nvect, err, 'k.', 'markersize', 24, 'linewidth', 2 )

%make plot pretty
title( 'Maximum error' ,'interpreter', 'latex','fontsize', 16)
xlabel( '$n$', 'interpreter', 'latex', 'fontsize', 16)
ylabel( '$||f(x) - f_a(x)||_2$', 'interpreter', 'latex', 'fontsize', 16)

set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )

set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [15 15])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 15 15])
set(gcf, 'PaperPosition', [0 0 15 15])

svnm = 'error';
print( '-dpng', svnm, '-r200' )

