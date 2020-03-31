%% AE 370 HW 1
%
% Date:     01/31/2020
% Modified: 02/01/2020
% Author:   Tiger Hou (linyih2)

%% Problem 1

%% Problem 2
close all hidden
clear;clc

% -----------------------------------------
% Part a

% -----------------------------------------
% Part b

% -------------------------
% Lagrange basis

n_vect  = [3,10,50,200];
f       = @(x) sin(cos(3*x)); 
intv    = [1,4];
err_L   = zeros(size(n_vect));

for i = 1:length(n_vect)
    n = n_vect(i);
    
    xp = (intv(1)+intv(2))/2 + (intv(1)-intv(2))/2*cos((0:n)*pi/n);
    fp = f(xp);
    dp = fp;
    A = eye(n+1);
    
    xx = linspace(intv(1),intv(2),1000);
    pn_L = 0;
%     L = @(x) 0;
    for j = 1:n+1
        Li = @(x) 1;
        
        for k = 1:n+1
            if k ~= j
                Li = @(x) (x-xp(k)) ./ (xp(j)-xp(k)) .* Li(x);
            end
        end
%         L = @(x) L(x) + dp(j) .* Li(x);
        pn_L = pn_L + dp(j) * Li(xx);
    end
    
    figure(1)
    subplot(4,1,i)
    hold on
    plot(xx,f(xx),'b-','lineWidth',2)
    plot(xp,f(xp),'k.','markersize',20)
    plot(xx,pn_L,'r--','linewidth',3)
    title([num2str(n+1) ' Sample Points'])
    xlabel('x')
    ylabel('f(x)')
    err_L(i) = norm(f(xx)-pn_L);
end

hold off

% -------------------------
% Monomial basis

err_M   = zeros(size(n_vect));

for i = 1:length(n_vect)
    n = n_vect(i);
    
    xp = (intv(1)+intv(2))/2 + (intv(1)-intv(2))/2*cos((0:n)*pi/n);
    fp = f(xp);
    
    A = zeros(n+1,n+1);
    for j = 1:n+1
        A(:,j) = xp' .^ (j-1);
    end
    dp = A \ (fp');
    
    xx = linspace(intv(1),intv(2),1000);
    M = @(x) 0;
    for j = 1:n+1
        M = @(x) dp(j) .* x.^(j-1) + M(x);
    end
    pn_M = M(xx);
    
    figure(1)
    hold on
    subplot(4,1,i)
    plot(xx,pn_M,'g-','lineWidth',1)
    legend('f(x)','sample points','fa(x), Lagrange','fa(x), Monomial')
    title([num2str(n+1) ' Sample Points'])
    xlabel('x')
    ylabel('f(x)')
    err_M(i) = norm(f(xx)-pn_M);
end

sgtitle('Lagrangian and Monomial Approximations of f(x) = sin(cos(3x))')

hold off

figure(100)
semilogy(n_vect,err_L,'r.','markersize',24,'linewidth',2), hold on
semilogy(n_vect,err_M,'g.','markersize',24,'linewidth',2), hold off
legend('Lagrangian Basis','Monomial Basis')
title('Cumulative Error of Approximation Functions')
xlabel('Sample Size')
ylabel('Error')


% -----------------------------------------
% Part c

n_vect  = [3,6,12,24,48,96];
intv    = [1,4];
condn_L = zeros(length(n_vect),1);
condn_M = condn_L;

% Lagrangian Basis
for i = 1:length(n_vect)
    n = n_vect(i);
    A = eye(n+1);
    condn_L(i) = cond(A);
end

% Monomial Basis
for i = 1:length(n_vect)
    n = n_vect(i);
    xp = (intv(1)+intv(2))/2 + (intv(1)-intv(2))/2*cos((0:n)*pi/n);
    A = zeros(n+1,n+1);
    for j = 1:n+1
        A(:,j) = xp' .^ (j-1);
    end
    condn_M(i) = cond(A);
end

figure(200)
semilogy(n_vect,condn_L,'r.','markersize',24,'linewidth',2), hold on
semilogy(n_vect,condn_M,'g.','markersize',24,'linewidth',2), hold off
legend('Lagrangian Basis','Monomial Basis')
title('Condition Number as a Function of Sample Size')
xlabel('n')
ylabel('Condition Number')

%% Problem 3
close all hidden
clear;clc

% -------------------------
% Lagrange basis

n_vect  = [5,10,15,20];
f       = @(x) 1 ./ (1+x.^2); 
intv    = [-5,5];

for i = 1:length(n_vect)
    n = n_vect(i);
    
    xp_E = intv(1):(intv(2)-intv(1))/n:intv(2);
    xp_C = (intv(1)+intv(2))/2 + (intv(1)-intv(2))/2*cos((0:n)*pi/n);
    fp_E = f(xp_E);
    fp_C = f(xp_C);
    dp_E = fp_E;
    dp_C = fp_C;
    
    xx = linspace(intv(1),intv(2),1000);
    pn_E = 0;
    for j = 1:n+1
        Li = @(x) 1;
        
        for k = 1:n+1
            if k ~= j
                Li = @(x) (x-xp_E(k)) ./ (xp_E(j)-xp_E(k)) .* Li(x);
            end
        end
        pn_E = pn_E + dp_E(j) * Li(xx);
    end
    
    pn_C = 0;
    for j = 1:n+1
        Li = @(x) 1;
        
        for k = 1:n+1
            if k ~= j
                Li = @(x) (x-xp_C(k)) ./ (xp_C(j)-xp_C(k)) .* Li(x);
            end
        end
        pn_C = pn_C + dp_C(j) * Li(xx);
    end
    
    figure(3)
    subplot(4,1,i)
    hold on
    plot(xx,f(xx),'b-','lineWidth',2)
    plot(xp_E,f(xp_E),'ro','markersize',7)
    plot(xp_C,f(xp_C),'go','markersize',7)
    plot(xx,pn_E,'r-','linewidth',2)
    plot(xx,pn_C,'g-','linewidth',2)
    legend( 'f(x)',...
            'equidistant points',...
            'Chebyshev points',...
            'fa(x), equidistant',...
            'fa(x), Chebyshev')
    title(['n = ' num2str(n)])
    xlabel('x')
    ylabel('f(x)')
    
end

sgtitle('Lagrange Approximation Using Equidistant and Chebyshev Points')

hold off