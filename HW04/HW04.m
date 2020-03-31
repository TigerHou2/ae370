%% Problem 1
close all hidden
clear;clc
format long

% a) Equispaced Points with Lagrange Basis

nvect = [5, 15, 30, 50];

fcn = @(x) 1 ./ (1 + x.^2);

%exact integral
xs = sym( 'x' );
%Can use int command to analytically integrate fcn over domain [-5,5]
int_exact = int( fcn( xs ), xs, -5, 5 );

err = zeros( length( nvect ), 1 );

for j = 1 : length( nvect )
    
    n = nvect( j );
    
    xj = linspace(-5,5,n+1);
    
    %define Lagrange basis vectors
    intval = 0;
    for i = 1 : n + 1
        L_i = 1;
        
        %vector indexing can't start at zero, so go from 1 to n+1
        for k = 1 : n + 1
            if k ~= i
                L_i = ( xs - xj( k ) )./( xj( i ) - xj( k ) ) .* L_i;
            end
        end
        
        %Again, use int (this time to compute integral of Lagrange fcn)
        Li_int = int( L_i, xs, -5, 5 );
        
        %Update integral
        intval = intval + Li_int * fcn(xj(i));
        
    end
    
    err( j ) = abs( int_exact - intval );
    
end

%plot error
figure(100)
semilogy( nvect, err, 'k.', 'markersize', 26 )
grid on


% b) Chebyshev Points with Lagrange Basis

err = zeros( length( nvect ), 1 );

for j = 1 : length( nvect )
    
    n = nvect( j );
    
    xj = 5*cos( (0:n)*pi/n );
    
    %define Lagrange basis vectors
    intval = 0;
    for i = 1 : n + 1
        %How to get Li???
        L_i = 1;
        
        %vector indexing can't start at zero, so go from 1 to n+1
        for k = 1 : n + 1
            if k ~= i
                L_i = ( xs - xj( k ) )./( xj( i ) - xj( k ) ) .* L_i;
            end
        end
        
        Li_int = int( L_i, xs, -5, 5 );
                
        intval = intval + Li_int * fcn(xj(i));
        
    end
    
    err( j ) = abs( intval - int_exact );
    
end

%plot error
figure(100), hold on
semilogy( nvect, err, 'r.', 'markersize', 26 )
grid on


% c) Equispaced Points with Trapezoid Rule

err = zeros( length( nvect ), 1 );

for j = 1 : length( nvect )
    
    n = nvect( j );
    
    xj = linspace(-5,5,n+1);
    
    dx = xj( 2 ) - xj( 1 ); %spacing between points
    
    %do end points first
    intval = dx / 2 * (fcn( xj( 1 ) ) + fcn( xj( n+1 ) ));
    
    %remaining points
    for jj = 2 : n
        
        intval = intval + dx * fcn( xj( jj ) );
        
    end
    
    err( j ) = abs( intval - int_exact );
    
end

%plot error
figure(100), hold on
semilogy( nvect, err, 'b.', 'markersize', 26 )
grid on


%make plot pretty
h = legend( 'Global, uniform points', 'Global, Cheb points', 'Comp. Trap' );
set( h, 'location', 'SouthWest', 'interpreter', 'latex', 'fontsize', 16)
xlabel( '$n$', 'interpreter', 'latex', 'fontsize', 16)
ylabel( 'error', 'interpreter', 'latex', 'fontsize', 16)

set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )

set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [15 15])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 15 15])
set(gcf, 'PaperPosition', [0 0 15 15])

svnm = 'error_compare';
print( '-dpng', svnm, '-r200' )


%% Problem 3

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
    
    xj = 0 : dx : 2;
    
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


%% Problem 4

close all hidden
clear;clc
format long

%Assumes a vector intvalh, of length length(nvect), is computed from 
%part a) using the Trapezoid rule with length h between points. Part b) 
%will compute the trapezoid rule at length h/2 and combine the two to 
%produce an improved estimate of the integral.

nvect = [4, 16, 24, 48];

fcn = @(x) (x.^2).*sin(10.*x);

%exact answer
xs = sym( 'x' );
int_exact = int( fcn( xs ), xs, 0, 2 );


% a) Equispaced Points with Trapezoid Rule

err = zeros( length( nvect ), 1 );
intvalh = zeros( length( nvect ), 1 );

for j = 1 : length( nvect )
    
    n = nvect( j );
    
    dx = 2/n; %spacing between points
    
    xj = 0 : dx : 2;
    
    %do end points first
    intval = dx / 2 * (fcn( xj( 1 ) ) + fcn( xj( n + 1 ) ));
    
    %remaining points
    for jj = 2 : n
        
        intval = intval + dx * fcn( xj( jj ) );
        
    end
        
    err( j ) = abs( int_exact - intval );
    intvalh( j ) = intval;
    
end

%plot error
figure(100)
semilogy( nvect, err, 'b.', 'markersize', 26 )


% b) Richardson Extrapolation by adding to Part a)

err = zeros( length( nvect ), 1 );

%--evaluate trap rule at h/2 and combine this with result from a) to get
%  Richardson extrapolated value.
intvalhb2 = zeros( length( nvect ), 1 );

for j = 1 : length( nvect )
    
    n = nvect( j );
    
    %h/2
    dx = 2/n / 2; %spacing between points
    
    %corresponding x points
    xj = 0 : dx : 2;
        
    %do end points first
    intvalhb2( j ) = dx / 2 * (fcn( xj( 1 ) ) + fcn( xj( 2*n + 1 ) ));
    
    %remaining points
    for jj = 2 : 2*n
        
        intvalhb2( j ) = intvalhb2( j ) + dx * fcn( xj( jj ) );
        
    end
    
    %Richardson extrapolated value
    int_rich = ( 4 * intvalhb2( j ) - intvalh( j ) ) / 3;
    
    err( j ) = abs( int_rich - int_exact );
    
end


%plot error
figure(100), hold on
semilogy( nvect, err, 'r.', 'markersize', 26 )
grid on

%make plot pretty
h = legend( 'Composite Trapezoid', 'Richardson Extrapolation' );
set( h, 'location', 'SouthWest', 'interpreter', 'latex', 'fontsize', 16)
xlabel( '$n$', 'interpreter', 'latex', 'fontsize', 16)
ylabel( 'max error', 'interpreter', 'latex', 'fontsize', 16)

set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )

set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [15 15])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 15 15])
set(gcf, 'PaperPosition', [0 0 15 15])

svnm = 'error_q4';
print( '-dpng', svnm, '-r200' )