%% Problem 1

close all
clear;clc

%solve -u'' + 2*u = f over 0 < x < 1
%with BCs u(0) = u(3) = 0;

%left and right bounds
xl = 0; xr = 3;
ln = xr - xl;

%exact sol
uex = @(x) sin( 2* sin(4*pi*(x-xl)/ln) );


%f(x) =
fcn = @(x) (32*pi^2*cos(2*sin((4*pi*(x-xl))/ln)).*sin((4*pi*(x-xl))/ln))/ln^2 ...
    +(64*pi^2*sin(2*sin((4*pi*(x-xl))/ln)).*cos((4*pi*(x-xl))/ln).^2)/ln^2 + ...
    2*(sin( 2* sin(4*pi*(x-xl)/ln) ));


%# of n points to use
nvect = [40; 80; 120; 240];

%initialize error vect
err = zeros( size( nvect ) );

for j = 1 : length(nvect)

    n = nvect( j ) ;
    dx = ( xr - xl ) / n;

    %define grid (including boundary pts)
    xx = linspace(xl,xr,n+1)';

    %initialize soln vector and rhs vector
    u = zeros( n-1,1 );
    b = u;

    %Build G (matrix associated with (10) from partial solutions)
    G = 1*dx/3 * ( 4*eye(n-1) + diag(ones(n-2,1),1) + diag(ones(n-2,1),-1) ) ...
        + 1/dx * ( 2*eye(n-1) - diag(ones(n-2,1),1) - diag(ones(n-2,1),-1) );

    %Build RHS vector, b, by looping through each element
    for jj = 1 : n-1


        %--RHS (be very careful with indexing here!)

            %define x from j-1 to j
            x_lft = (xl+(jj-1)*dx):dx:(xl+(jj)*dx);
            %define phij from j-1 to j (upward sloping part of phij)
            phij_lft = 1/dx * ( x_lft-xl-(jj-1)*dx );
            %inner product contribution from j-1 to j
                %can use Matlab's trapz for simplicity, but other quadrature
                %rules are OK too!
            b(jj) = trapz( x_lft, fcn( x_lft ).* phij_lft );

            %define x from j to j+1 (downward sloping part of phij)
            x_rgt = (xl+(jj)*dx):dx:(xl+(jj+1)*dx);
            %phij from j to j+1
            phij_rgt = -1/dx * ( x_rgt-xl-(jj+1)*dx );
            %inner product contribution from j-1 to j
            b(jj) = b(jj) + trapz( x_rgt, fcn( x_rgt ).* phij_rgt );

        %--

    end

    %Solve for uhat
    uhat = G\b;

    %Add BCs to ua
    uhat = [uex(xl); uhat; uex(xr)];

    %--plot soln
        figure(1),subplot(2,2,j)
        plot( xx, uex(xx), 'k-', 'linewidth', 2 ), hold on
        plot( xx, uhat, 'r--', 'linewidth', 2 ), hold on

        %make plot pretty
        title( ['$n = ', num2str( n ),'$'] ,'interpreter', 'latex',...
            'fontsize', 16)
        if j == 1
            h = legend( '$u(x)$', '$\hat{u}(x)$');
        end

        if j <= 2
            set( gca, 'XTick', [] )
        else
            xlabel( '$x$', 'interpreter', 'latex', 'fontsize', 16)

        end
        set(h, 'location', 'NorthEast', 'Interpreter', 'Latex', 'fontsize', 16 )
        set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )
        axis([0 3 -1.5 1.5])

        set(gcf, 'PaperPositionMode', 'manual')
        set(gcf, 'Color', [1 1 1])
        set(gca, 'Color', [1 1 1])
        set(gcf, 'PaperUnits', 'centimeters')
        set(gcf, 'PaperSize', [25 25])
        set(gcf, 'Units', 'centimeters' )
        set(gcf, 'Position', [0 0 25 25])
        set(gcf, 'PaperPosition', [0 0 25 25])
    %--

    %error
    err( j ) = norm( uhat - uex(xx) ) / norm( uex(xx) );


end

svnm = 'q1_plot';
print( '-dpng', svnm, '-r200' )

%plot err
figure(100)

%plot dx^2 line to show err scales correctly
c = err(end)/(dx^2);
loglog( nvect, c*(ln./nvect).^2, 'k--', 'linewidth', 2 ), hold on

%plot err
loglog( nvect, err , 'b.', 'markersize', 26 )
xlim([10 1000])

%make pretty
xlabel( '$n$', 'interpreter', 'latex', 'fontsize', 16)
ylabel( '$||\textbf{e}||/||\textbf{u}_e||$ ', 'interpreter', 'latex', 'fontsize', 16)

set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )

set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [15 15])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 15 15])
set(gcf, 'PaperPosition', [0 0 15 15])

svnm = 'q1_error';
print( '-dpng', svnm, '-r200' )
