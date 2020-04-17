clear all, close all


%%

%solve -u'' + 2*u = f over 0 < x < 1
%with BCs u(0) = u(3) = 0;

%left and right bounds
xl = 0; xr = 3;
ln = xr - xl;

%exact sol
uex = @(x) sin( 2* sin(4*pi*(x-xl)/ln) );


%f(x) = 
fcn = @(x)  (32*pi^2*cos(2*sin((4*pi*(x-xl))/ln)).*sin((4*pi*(x-xl))/ln))/ln^2 ...
    +(64*pi^2*sin(2*sin((4*pi*(x-xl))/ln)).*cos((4*pi*(x-xl))/ln).^2)/ln^2 + ...
    2*(sin( 2* sin(4*pi*(x-xl)/ln) ));


%# of n points to use
nvect = [10; 20; 40; 60];

%initialize error vect
err = zeros( size( nvect ) );

%x points on which to evaluate u
xx = linspace(0,3,1000)';

for j = 1 : length( nvect )
    
    n = nvect(j);
    
    xj = linspace(xl,xr,n)';
    
    %Store coefficients in c vector
    c = zeros( n,1);
    for jj = 1 : n
        sin_basis_jj = @(x) sin(jj*pi*(x-xl)/ln);
        
        c(jj) = ip_s(fcn,sin_basis_jj,xj)/ip_E(sin_basis_jj,sin_basis_jj,xj);
    end
        

    %Represent u on the finer xx domain
    u = 0;
    for kk = 1 : n
        sin_basis_kk = @(x) sin(kk*pi*(x-xl)/ln);
        
%         u = u + c(kk) * sin_basis_kk(xx);
        
        u = u + c(kk) * sin((kk*pi*(xx-xl))/(ln));

    end
    
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
    err( j ) = norm( u - uex(xx) )/norm(uex(xx)) ;
    
    
end

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





function ip_eval = ip_s(func,basis,x_points)
    integrand = @(x) func(x).*basis(x);
    
    integrand_eval = integrand(x_points);
    
    h = x_points(2) - x_points(1);
    n = length(x_points);
    
    %do end points first
    ip_eval = (h/2)*(integrand_eval(1) + integrand_eval(end));
    
    %remaining points
    for jj = 2 : (n-1)
        
        ip_eval = ip_eval + h*integrand_eval(jj);
        
    end
    
    
%     ip_eval = trapz(integrand(x_points),x_points);
%     ip_eval = integral(integrand,x_points(1),x_points(end));
end

function ip_eval = ip_E(basis1,basis2,x_points)
    syms x
    
    basis1_d1 = matlabFunction(diff(basis1,x));
    basis2_d1 = matlabFunction(diff(basis2,x));
    
    integrand = @(x) basis1_d1(x).*basis2_d1(x) + 2.*basis1(x).*basis2(x);
    
    integrand_eval = integrand(x_points);
    
    
    h = x_points(2) - x_points(1);
    n = length(x_points);
    
    %do end points first
    ip_eval = (h/2)*(integrand_eval(1) + integrand_eval(end));
    
    %remaining points
    for jj = 2 : (n-1)
        
        ip_eval = ip_eval + h*integrand_eval(jj);
        
    end
    
%     ip_eval = trapz(integrand(x_points),x_points);
%     ip_eval = integral(integrand,x_points(1),x_points(end));
end


