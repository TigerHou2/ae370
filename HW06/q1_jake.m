clear all, close all, clc


%%

%solve u'' + 3*u' = f over 0 < x < 1
%with BCs u(0) = u(1) = 0;

%use central differences

%exact sol
uex = @(x) sin( 5.*pi*x ).*sin( 3*pi.*x );

%f(x) = 
fcn = @(x) 30*pi^2.*cos(3.*pi.*x).*cos(5.*pi.*x) - ...
    34*pi^2*sin(3.*pi.*x).*sin(5.*pi.*x) + ...
    3*(3.*pi.*cos(3.*pi.*x).*sin(5.*pi.*x) + ...
    5.*pi.*cos(5.*pi.*x).*sin(3.*pi.*x) );

%left and right bounds
xbds = [0 1];

%# of n points to use
nvect = [20,40,80,160];

%initialize error vect
err = zeros( size( nvect ) );

for j = 1 : length( nvect )
    
    n = nvect( j) ;
    dx = (xbds(2) - xbds(1))/n;
    
    %x points
    xx = xbds(1) + dx*((1:(n+1)) - 1)';
    
    %--LHS matrices (2nd deriv &  1st deriv matrix)
    
        %**second deriv
            %interior points

            %Central Difference
            L = zeros(n+1);
            for i = 2:n
                L(i,i-1:i+1) = [1/dx^2 , -2/dx^2 , 1/dx^2];
            end

        %**
        
        %**first deriv
            %interior points

            %Central Difference
            D = zeros(n+1);
            for k = 2:n
                D(k,k-1:k+1) = [-1/(2*dx) , 0 , 1/(2*dx)];
            end

        %**
        
        %**combine and incorporate BCs
             
            A = L + 3*D;

            %BCs
            %remove erroneous entries in first & last rows
            A(1,:) = 0; 
            A(n+1, :) = 0; 
            
            %correct entries in first & last rows
            A(1,1) = 1;
            A(n+1,n+1) = 1;
        %**
        
    %--
    
    %rhs vector
    f = [0 ; fcn(xx(2:(end-1))) ; 0];
    
    %approximate solution
    u = A \ f;
    
    %--plot soln
        figure(1),subplot(2,2,j)
        plot( xx, uex(xx), 'k-', 'linewidth', 2 ), hold on
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
    err( j ) = sqrt(dx)*norm( u - uex(xx) );
    
    
end

% svnm = 'n_plots_q1';
% print( '-dpng', svnm, '-r200' )


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

set(gca, 'TickLabelInterpreter','latex', 'fontsize', 16 )

set(gcf, 'PaperPositionMode', 'manual')
set(gcf, 'Color', [1 1 1])
set(gca, 'Color', [1 1 1])
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperSize', [15 15])
set(gcf, 'Units', 'centimeters' )
set(gcf, 'Position', [0 0 15 15])
set(gcf, 'PaperPosition', [0 0 15 15])


% svnm = 'error_q1';
% print( '-dpng', svnm, '-r200' )

