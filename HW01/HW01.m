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

n_vect = [3,10,50,200];
f = @(x) sin(cos(3*x)); 
intv = [1,4];
err = zeros(size(n_vect));
condn = zeros(size(n_vect));

for i = 1:length(n_vect)
    n = n_vect(i);
    
    xp = (intv(1)+intv(2))/2 + (intv(1)-intv(2))/2*cos((0:n)*pi/n);
    fp = f(xp);
    dp = fp;
    A = eye(n+1);
    
    xx = linspace(intv(1),intv(2),1000);
    pn = 0;
    for j = 1:n+1
        Li = @(x) 1;
        
        for k = 1:n+1
            if k ~= j
                Li = @(x) (x-xp(k)) ./ (xp(j)-xp(k)) .* Li(x);
            end
        end
        pn = pn + dp(j) * Li(xx);
    end
    
    figure(1)
    subplot(2,2,i)
    plot(xx,f(xx),'b-','lineWidth',2), hold on
    plot(xp,f(xp),'k.','markersize',20)
    plot(xx,pn,'r--','linewidth',3)
    
    err(i) = norm(f(xx)-pn);
    condn(i) = cond(A);
end

hold off

figure(101)
semilogy(n_vect,err,'k.','markersize',24,'linewidth',2)
figure(201)
semilogy(n_vect,condn,'k.','markersize',24,'linewidth',2)


% -------------------------
%% Monomial basis
close all
clear;clc

n_vect = [3,10,50,200];
f = @(x) sin(cos(3*x)); 
intv = [1,4];
err = zeros(size(n_vect));
condn = zeros(size(n_vect));

for i = 1:length(n_vect)
    n = n_vect(i);
    
    xp = (intv(1)+intv(2))/2 + (intv(1)-intv(2))/2*cos((0:n)*pi/n);
    fp = f(xp);
    
    A = zeros(n+1,n+1);
    for j = 1:n+1
        A(:,j) = xp' .^ (j-1);
    end
    dp = A\(fp');
    
    xx = linspace(intv(1),intv(2),1000);
    pn = 0;
    Li = @(x) 0;
    for k = 1:n+1
        Li = @(x) dp(k) .* x.^(k-1) + Li(x);
    end
    pn = Li(xx);
    
    figure(2)
    subplot(2,2,i)
    plot(xx,f(xx),'b-','lineWidth',2), hold on
    plot(xx,pn,'r--','linewidth',2)
    plot(xp,f(xp),'k.','markersize',24)
    
    err(i) = norm(f(xx)-pn);
    condn(i) = cond(A);
end

hold off

figure(102)
semilogy(n_vect,err,'k.','markersize',24,'linewidth',2)
figure(202)
semilogy(n_vect,condn,'k.','markersize',24,'linewidth',2)

% -----------------------------------------
% Part c
