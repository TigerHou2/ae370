close all hidden
clear;clc
format long

nvect = [4, 16, 24, 48];

fcn = @(x) sin(10.*pi*x);

%exact answer
xs = sym( 'x' );
int_exact = int( fcn( xs ), xs, 0, 2 );

% a) Equispaced Points with Trapezoid Rule

err = zeros( length( nvect ), 1 );
% intvalh = zeros( length( nvect ), 1 );

for j = 1 : length( nvect )
    
    n = nvect( j );
    
    dx = 2/n; %spacing between points
    
    xj = linspace(0,2,n+1);
    
    %do end points first
    intval = dx / 2 * (fcn( xj( 1 ) ) + fcn( xj(n+1) ));
    
    %remaining points
    for jj = 2 : n
        
        intval = intval + dx * fcn( xj( jj ) );
        
    end
        
    err( j ) = abs( int_exact - intval );
    
end

%plot error
figure(100)
semilogy( nvect, err, 'b.', 'markersize', 26 )
grid on
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )
% ylim([10^(-20) 1])
xlabel( '$n$', 'interpreter', 'latex', 'fontsize', 16)
ylabel( 'error', 'interpreter', 'latex', 'fontsize', 16)

set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [15 15])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 15 15])
set(gcf, 'PaperPosition', [0 0 15 15])

svnm = 'error_q3';
print( '-dpng', svnm, '-r200' )