%% main.m
%
% Author: 
%   Tiger Hou
%
% Description:
%   Finite element method for reflected waves in changing physical media
%
%% simulation
%
% https://www-users.math.umn.edu/~olver/num_/lnp.pdf
%
close all
clear;clc

%solve the wave eqn over 0 < x < 1
%with BCs u(10,t) = u(20,t) = 0
%and IC u(x,0) = exp(-10*x);

saveGIF = false;

%params for problem
a = 0; 
b = 1;
c = @(x) 0.25*(2-atan(40*(x-b/2)));
ln = b-a;

%initial condition
u0 = @(x) exp(-200.*(x-0.3).^2);
v0 = @(x) 0;

%# of n points to use
nvect = 400;

T = 4; % final time
dt = (b-a)/nvect(end)/5;

%initialize u vect
uvect = cell(size(nvect));

% iterate through interior point sizes
for j = 1 : length( nvect )

    % Build n, xj points, A matrix and g vector
        % # of grid points
        n = nvect( j );

        % build interp points
        xj = linspace(a,b,n+1)';

        %grid spacing (uniform)
        dx = ( b - a ) / n;
        
        %  sigma function
        s = c(xj) * dt / dx;
        
        % build A matrix
        diag_ct = diag(2*(1-s(2:end-1).^2));
        diag_dn = diag(s(2:end-2).^2,-1);
        diag_up = diag(s(3:end-1).^2, 1);
        A = diag_ct + diag_dn + diag_up;

        % also build identity mat (same size as A)
        I = eye(size(A));

        % build b for this set of xj
        g = @(t) 0 * xj(2:end-1) * t;

    % Iterate through time
    
        % build time vector
        tvect = 2*dt:dt:T;
        
        % second order states initialization
        um1 = u0(xj(2:end-1));
        uu  = 0.5 * A * um1 + dt * v0(xj(2:end-1)) + 0.5 * g(0);
        uuvect = repmat(uu,1,length(tvect));
        uuvect(:,1) = um1;
        uuvect(:,2) = uu;
        
        for i = 3:length(tvect)
            
            up1 = A * uu - um1 + g(tvect(i));
            um1 = uu;
            uu  = up1;
            
            uuvect(:,i) = up1;
            
        end
        
        uvect{j} = uuvect;

end


%% plotting

temp = uvect{end};
xx = linspace(a,b,size(temp,1));
wave = line(NaN,NaN,'LineWidth',1,'Color','k');

figure(1)

hold on

% make background gradient
yy = linspace(min(temp(:)),max(temp(:)),100);
[X,Y] = meshgrid(xx,yy);
Z = c(X);
Z = Z / max(Z(:));
surf(X,Y,Z,'EdgeColor','none'...
          ,'FaceColor','interp'...
          ,'FaceAlpha',0.2)
cmap = cool(100);
colormap(cmap(ceil(min(Z(:))*100):100,:))
view(2)

% animate and make gif
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'test.gif';
lim_x = [a;b];
lim_y = [min(temp(:));max(temp(:))];
lim_z = [0;1];
update_view(lim_x,lim_y,lim_z);
playtime = T;
pausetime = playtime/size(temp,2);
frameRate = 1/pausetime;

for i = 1:size(temp,2)

    set(wave,'XData',xx,'YData',temp(:,i))
    
    if saveGIF
        
        % Capture the plot as an image 
        frame = getframe; 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        % Write to the GIF File 
        if i == 1 
          imwrite(imind,cm,filename,'gif','Loopcount',inf,...
                                          'DelayTime',pausetime); 
        else 
          imwrite(imind,cm,filename,'gif','WriteMode','append',...
                                          'DelayTime',pausetime); 
        end 
    end
    
    pause(pausetime)
    
end
hold off

%% supporting functions

function update_view(lim_x,lim_y,lim_z)

    xlim(lim_x)
    ylim(lim_y)
    zlim(lim_z)

end