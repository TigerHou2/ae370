close all
clear;clc

syms a ii x dx real

%% phi(i-1), phi(i)

phi_m1 = -1/dx*(x-a-ii*dx);
phi =     1/dx*(x-a-(ii-1)*dx);

lb = a+(ii-1)*dx;
ub = a+ii*dx;

IP =     int( diff(phi_m1,x)*diff(phi,x), x, lb, ub ) ...
   + 2 * int( phi_m1*phi                , x, lb, ub )

%% phi(i), phi(i)

phi_l =  1/dx*(x-a-(ii-1)*dx);
phi_u = -1/dx*(x-a-(ii+1)*dx);

lb = a+(ii-1)*dx;
mb = a+ii*dx;
ub = a+(ii+1)*dx;

IP =     int( diff(phi_l,x)*diff(phi_l,x), x, lb, mb ) ...
   + 2 * int( phi_l*phi_l                , x, lb, mb ) ...
   +     int( diff(phi_u,x)*diff(phi_u,x), x, mb, ub ) ...
   + 2 * int( phi_u*phi_u                , x, mb, ub )

%% phi(i), phi(i+1)

phi =    -1/dx*(x-a-(ii+1)*dx);
phi_p1 =  1/dx*(x-a-ii*dx);

lb = a+ii*dx;
ub = a+(ii+1)*dx;

IP =     int( diff(phi,x)*diff(phi_p1,x), x, lb, ub ) ...
   + 2 * int( phi*phi_p1                , x, lb, ub )