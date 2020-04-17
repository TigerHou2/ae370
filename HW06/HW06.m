%% Problem 1

close all
clear;clc

%solve u'' + 3*u' = f over 0 < x < 1
%with BCs u(0) = u(1) = 0;

%use central differences

%exact sol
uex = @(x) sin( 5.*pi*x ).*sin( 3*pi.*x );

%f(x) = 
fcn = @(x) 30 * pi^2.*cos(3.*pi.*x).*cos(5.*pi.*x) ...
         - 34 * pi^2.*sin(3.*pi.*x).*sin(5.*pi.*x) ...
         + 3  * ( 3.*pi.*cos(3.*pi.*x).*sin(5.*pi.*x) ...
                 +5.*pi.*cos(5.*pi.*x).*sin(3.*pi.*x) );

%left and right bounds
xbds = [0 1];

%# of n points to use
nvect = [20, 40, 80, 160];

%initialize error vect
err = zeros( size( nvect ) );

for j = 1 : length( nvect )
    
    n = nvect( j ) ;
    dx = ( xbds(2) - xbds(1) ) / n;
    
    %x points
    xx = linspace(xbds(1),xbds(2),n+1)';
    
    %--LHS matrices (2nd deriv &  1st deriv matrix)
    
        %**second deriv. interior points (Central Difference)
        L = (1/dx^2) * (diag(ones(n,1),1) + diag(ones(n,1),-1) - 2*eye(n+1));

        %** first deriv. interior points (Central Difference)
        D = (1/dx) * (0.5 * diag(ones(n,1),1) - 0.5 * diag(ones(n,1),-1));

        %**combine and incorporate BCs
        A = L + 3*D;
    
    %--

    %remove erroneous entries in first & last rows
    A(1,:) = 0; 
    A(n+1,:) = 0; 

    %correct entries in first & last rows
    A(1,1) = 1;
    A(end) = 1;
    
    %rhs vector
    f = fcn(xx);
    f(1) = 0;
    f(end) = 0;
    
    %approximate solution
    u = A \ f;
    
    % exact solution
    ue = uex(xx);
    
    %--plot soln
        figure(1),subplot(2,2,j)
        plot( xx, ue, 'k-', 'linewidth', 2 ), hold on
        plot( xx, u , 'r--', 'linewidth', 2 ), hold on

        %make plot pretty
        title( ['$n = ', num2str( n ),'$'] ,'interpreter', 'latex',...
            'fontsize', 16)
        if j == 1
            h = legend( '$\overline{u}(x)$', '$u(x)$');
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
       
    %--

    %error
    err( j ) = sqrt(dx)*norm( u - ue );
    
    
end

svnm = 'q1_plot';
print( '-dpng', svnm, '-r200' )

%plot err
figure(100)

%plot dx^2 line to show err scales correctly
c = err(end)/(dx^2);
loglog( nvect, c./((nvect+1).^2), 'k--', 'linewidth', 2 ), hold on

%plot err
loglog( nvect, err , 'b.', 'markersize', 26 )
xlim([10 500])

h = legend( '$O(\Delta x^2)$', '$||e||_2$');
set(h, 'location', 'NorthEast', 'Interpreter', 'Latex', 'fontsize', 16 )

%make pretty
xlabel( '$n$', 'interpreter', 'latex', 'fontsize', 16)
ylabel( '$||\textbf{e}||_2$ ', 'interpreter', 'latex', 'fontsize', 16)

grid(gca,'minor')
grid on

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


%% Problem 2

close all
clear;clc

%solve -u'' + 2*u = f over 0 < x < 3
%with BCs u(0) = u(3) = 0;

%left and right bounds
xl = 0; xr = 3;
ln = xr - xl;

%exact sol
uex = @(x) sin( 2* sin(4*pi*(x-xl)/ln) );


%f(x)
fcn = @(x) ...
    (1/ln^2) ...
        * ( (32*pi^2*cos(2*sin((4*pi*(x-xl))/ln)).*sin((4*pi*(x-xl))/ln)) ...
          + (64*pi^2*sin(2*sin((4*pi*(x-xl))/ln)).*cos((4*pi*(x-xl))/ln).^2) ) ...
  + 2*(sin( 2* sin(4*pi*(x-xl)/ln) ));

syms x
dfcn = matlabFunction(diff(fcn(x)));

%# of n points to use
nvect = [10; 20; 40; 60];

%initialize error vect
err = zeros( size( nvect ) );

%x points on which to evaluate u
xx = linspace(xl,xr,1000);

for j = 1 : length( nvect )
    
    n = nvect( j );
    
    %Store coefficients in c vector
    c = zeros(n,1);
    for jj = 1 : n
        
        c(jj) = trapz(xx,fcn(xx).*sin(jj.*pi.*(xx-xl)/ln)) / ...
          ( 1*trapz(xx, (cos(jj.*pi.*(xx-xl)/ln)*jj*pi/ln) .* (cos(jj.*pi.*(xx-xl)/ln)*jj*pi/ln) ) ...
           +2*trapz(xx, (sin(jj.*pi.*(xx-xl)/ln)) .* (sin(jj.*pi.*(xx-xl)/ln)) ) );

    end
    

    %Represent u on the finer xx domain
    u = 0;
    for kk = 1 : n
        
        u = u + c(kk) .* sin(kk.*pi.*(xx-xl)/ln);

    end
    
    ue = uex(xx);
    
    %--plot soln
        figure(1),subplot(2,2,j)
        plot( xx, ue, 'k-', 'linewidth', 2 ), hold on
        plot( xx, u, 'r--', 'linewidth', 2 ), hold on

        %make plot pretty
        title( ['$n = ', num2str( n ),'$'] ,'interpreter', 'latex',...
            'fontsize', 16)
        if j == 1
            h = legend( '$\overline{u}(x)$', '$u(x)$');
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
    %--

    %error
    err( j ) = norm( u - uex(xx) )/norm(uex(xx)) ;
    
    
end

svnm = 'q2_plot';
print( '-dpng', svnm, '-r200' )

%plot err
figure(100)

% %plot dx^2 line to show err scales correctly
% c = err(end)/(dx^2);
% loglog( nvect, c./((nvect+1).^2), 'k--', 'linewidth', 2 ), hold on

%plot err
loglog( nvect, err , 'b.', 'markersize', 26 )
xlim([10 500])

%make pretty
xlabel( '$n$', 'interpreter', 'latex', 'fontsize', 16)
ylabel( '$||\textbf{e}||_2$ ', 'interpreter', 'latex', 'fontsize', 16)

set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )

set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [15 15])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 15 15])
set(gcf, 'PaperPosition', [0 0 15 15])

svnm = 'q2_error';
print( '-dpng', svnm, '-r200' )