clear all, close all, clc

%solve heat eqn over 2 < x < 15
%with BCs u(2,t) = u(15,t) = 0
%and IC u(x,0) = 0
%This code does a convergence test in time

%params for problem
a = 2; 
b = 15;
kappa = 2; 
ln = b-a;
T = 5; %Final time to run to

%Make dx small in spatial convergence test so that spatial error
%doesn't pollute convergence
n = 3000; 
dx = (b-a)/n;  

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
dtvect = [2.5; 1; 5e-1; 2.5e-1; 1e-1];

%initialize error vect
err = zeros( size( dtvect ) );

for j = 1 : length( dtvect )


    %---Build xj points, A matrix and g vector
        %time step size
        dt = ????

        %Build interp points
        xj = ????

        %Build A matrix
        %Use truncated version from lecture notes
        A = ????

        %Also build identity mat (same size as A)
        I = ????

        %Build g for this set of xj
        g = @(t) ????

        %Build RHS for IVP, f(u,t)
        f = @(u,t) ?????
    %---

    %---Initialize for time stepping
        uk = ??????
        tk = 0;
        tvect = dt : dt : T;
    %---
    
    %---Do time stepping
    for jj = 1 : length( tvect )

        tkp1 = tk + dt;
        
        %Update solution at next time using trap method
        ukp1 = ????

        uk = ????
        tk = ????

    end
    %---
    
    err(j) = ?????
end


%--Error plots
    figure(2)
    c = err(end)*(1./dt^2);
    loglog( dtvect, c*(dtvect).^2, 'k--', 'linewidth', 2 ), hold on

    %plot err
    loglog( dtvect, err , 'b.', 'markersize', 26 )
    xlim([1e-2 100])
    
    
    %make pretty
    h = legend('$O(\Delta t^2)$', '$Error$');
    set(h, 'Interpreter','latex', 'fontsize', 16, 'Location', 'NorthWest' )
    
    xlabel( '$\Delta t$', 'interpreter', 'latex', 'fontsize', 16)
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




