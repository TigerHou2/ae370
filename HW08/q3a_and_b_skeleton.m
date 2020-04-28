clear all, close all, clc

%solve heat eqn over 2 < x < 15
%with BCs u(2,t) = u(15,t) = 0
%and IC u(x,0) = 0
%This code does a convergence test in space

%params for problem
a = 2; 
b = 15;
kappa = 2; 
ln = b-a;
T = 5; %Final time to run to
dt = 1e-3; %Make dt small in spatial convergence test so that time error 
           %doesn't pollute convergence


%exact sol
uex = @(x,t) sin( t.* sin(6.*pi.*(x-a)./ln) );

%initial condition
eta = @(x) uex(x,0);

%g(x) = 
fcn = @(x,t) - kappa*((36*t.^2.*pi^2.*cos((6.*pi.*(a - x))./ln).^2 ...
    .*sin(t.*sin((6.*pi.*(a - x))./ln)))./ln.^2 + ...
    (36.*t.*pi^2.*sin((6.*pi.*(a - x))./ln).*cos(t.*sin((6.*pi.*(a - x))...
    ./ln)))/ln.^2) - sin((6.*pi.*(a - x))./ln)...
    .*cos(t.*sin((6.*pi.*(a - x))./ln));

%# of n points to use
nvect = [20; 40; 80; 120; 240];

%initialize error vect
err = zeros( ???? );

for j = 1 : length( nvect )


    %---Build n, xj points, A matrix and g vector
        %# of grid points
        n = ????

        %Build interp points
        xj = ????

        %grid spacing (uniform)
        dx = ???
        
        %Build A matrix
        %Use truncated version from lecture notes
        A = ????

        %Also build identity mat (same size as A)
        I = eye(size(A));

        %Build g for this set of xj
        g = @(t) ???

        %Build RHS for IVP, f(u,t)
        f = @(u,t) ????
    %---

    %---Initialize for time stepping
        uk = eta(xj(2:end-1));
        tk = 0;
        tvect = dt : dt : T;

        %# snapshots to save (don't mess with this; it sets things up so
        %the solution is only saved a relatively small number of times 
        %to keep your data storage from growing to large)
        nsnps = 100;
        ind = max( 1, round(length(tvect)/nsnps) );
        tsv = tvect( 1 : ind : end );

        u = zeros( n-1, length(tsv));
        cnt = 1;
    %---
    
    %---Do time stepping
    for jj = 1 : length( tvect )

        tkp1 = tk + dt;
        
        %Update solution at next time using trap method
        ukp1 = ???

        %Update solution variable & time
        uk = ???
        tk = ???

        %Again, leave this. It sets things up to only save for a relatively
        %small # of times.
        if min(abs( tkp1-tsv ) ) < 1e-8

            u(:,cnt) = uk;
            cnt = cnt + 1;
        end

    end
    %---
    
    err(j) = norm( ukp1 - ???? ) / norm( uex(????) );

end


%--Waterfall plot of solns
    [X,T] = meshgrid( xj(2:end-1), tsv );

    figure(1), subplot(1,2,1)
    waterfall( X,T, u.' ), hold on
    set( gca, 'fontsize', 15, 'ticklabelinterpreter', 'latex' )
    title('$u_{FD}$', 'fontsize', 20, 'interpreter', 'latex')
    xlabel('$x$', 'fontsize', 15, 'interpreter', 'latex')
    ylabel('$t$', 'fontsize', 15, 'interpreter', 'latex')
    zlim([-1 1])

    subplot(1,2,2)
    waterfall( X,T, uex(X,T) ), hold on
    set( gca, 'fontsize', 15, 'ticklabelinterpreter', 'latex' )
    title('$u_{exact}$', 'fontsize', 20, 'interpreter', 'latex')
    xlabel('$x$', 'fontsize', 15, 'interpreter', 'latex')
    ylabel('$t$', 'fontsize', 15, 'interpreter', 'latex')
    zlim([-1 1])
    
    set(gcf, 'PaperPositionMode', 'manual')
    set(gcf, 'Color', [1 1 1])
    set(gca, 'Color', [1 1 1])
    set(gcf, 'PaperUnits', 'centimeters')
    set(gcf, 'PaperSize', [25 25])
    set(gcf, 'Units', 'centimeters' )
    set(gcf, 'Position', [0 0 25 25])
    set(gcf, 'PaperPosition', [0 0 25 25])
%--

%--Error plots
    figure(2)
    c = err(end)/(dx^2);
    loglog( ln./nvect, c*(ln./nvect).^2, 'k--', 'linewidth', 2 ), hold on

    %plot err
    loglog( ln./nvect, err , 'b.', 'markersize', 26 )
    xlim([1e-2 1])
    h = legend('$O(\Delta x^2)$', '$Error$');
    set(h, 'Interpreter','latex', 'fontsize', 16, 'Location', 'NorthWest' )
    
    %make pretty
    xlabel( '$\Delta x$', 'interpreter', 'latex', 'fontsize', 16)
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

%--




