%% main.m
%
% Author: 
%   Tiger Hou
%
% Description:
%   Finite difference method for reflected waves in changing physical media
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
play = true;

% define left/right boundaries
a = 0;     % left boundary
b = 3*pi;  % right boundary
ln = b-a;

% wave speed as a function of position
% --- switching function
c = @(x) (x< pi/2) * 1 + ...
         (x>=pi/2) * 1;
% --- sigmoid function
c = @(x) (3*pi/2-atan(15*(x-pi/2)))/pi;

% initial conditions
u0 = @(x) zeros(size(x));
v0 = @(x) zeros(size(x));

% boundary conditions
f = @(t) (-cos(t*2*pi)+1)*(pi/2-atan(10000*(t-1)))/pi;
g = @(t) 0;

% # of n points to use
nvect = [50, 100, 200, 500, 1000];

T = 8; % final time
dtvect = (b-a)./nvect/5; % dt must match dx in order of magnitude

% initialize u vect
uvect = cell(size(nvect));

% initialize error vector
final_vect = cell(size(nvect));

% iterate through interior point sizes
for j = 1 : length( nvect )

    % Build n, xj points, A matrix and g vector
        % # of grid points
        n = nvect( j );
        
        % choose time step matching dx
        dt = dtvect(j);

        % build interp points
        xj = linspace(a,b,n+1)';

        % grid spacing (uniform)
        dx = ( b - a ) / n;
        
        % sigma function
        s = c(xj) * dt / dx;
        
        % build A matrix
        diag_ct = diag(2*(1-s(2:end-1).^2));
        diag_dn = diag(s(2:end-2).^2,-1);
        diag_up = diag(s(3:end-1).^2, 1);
        A = diag_ct + diag_dn + diag_up;
        
        % special sigma functions for f and g
        ssf = c(a) * dt / dx;
        ssg = c(b) * dt / dx;

        % build h for this set of xj
        h = @(t) [ssf^2*f(t); zeros(size(xj,1)-4,1); ssg^2*g(t)];

    % Iterate through time
    
        % build time vector
        tvect = 2*dt:dt:T;
        
        % second order states initialization
        um1 = u0(xj(2:end-1));
        uu  = 0.5 * A * um1 + dt * v0(xj(2:end-1)) + 0.5 * h(0);
        uuvect = repmat(uu,1,length(tvect));
        uuvect(:,1) = um1;
        uuvect(:,2) = uu;
        
        for i = 3:length(tvect)
            
            up1 = A * uu - um1 + h(tvect(i));
            um1 = uu;
            uu  = up1;
            
            uuvect(:,i) = up1;
            
        end
        
        % store entire wave solution
        uvect{j} = uuvect;
        
        % store final state to compare errors later
        final_vect{j} = up1;

end

% --- error comparison
err = nan(length(nvect)-1,1);
xj_ref = xj;
uu_ref = up1;
for j = 1 : length( nvect ) - 1
    
    n = nvect(j);
    
    xj = linspace(a,b,n+1);
        
    [~,idx] = min( abs( xj_ref - xj(2:end-1) ) );
    
    uu = final_vect{j};
    
    err(j) = norm( uu_ref(idx) - uu );
    
end


%% plot error

figure(1)

cc = err(end)./(dx.^2);
loglog( ln./nvect(1:end-1), cc.*(ln./nvect(1:end-1)).^2, ...
        'k--', 'linewidth', 2 )
hold on

loglog( ln./nvect(1:end-1), err , 'b.', 'markersize', 26 )
legend('$O(\Delta x^2)$', 'Error', 'Location', 'Best');
hold off

grid(gca,'minor')
grid on


%% plot waves

figure(2)

% extract the highest resolution data
waveData = uvect{end};
% create x-axis, mapped to highest resolution data
xx = linspace(a,b,size(waveData,1));

% initialize wave
wave = line(NaN,NaN,'LineWidth',1,'Color','k');
hold on

% make background gradient
bgAlpha = 0.5;
yy = linspace(min(waveData(:)),max(waveData(:)),100);
X = meshgrid(xx,yy);
Z = c(X);
cmap = spring(100);
colormap(cmap(30:60,:))
bg = image(xx,yy,Z,'CDataMapping','scaled');
bg.AlphaData = bgAlpha;
bg_bar = colorbar;
bg_bar.TickLabelInterpreter = 'latex';
bg_bar.Label.String = 'Speed';
bg_bar.Label.Interpreter = 'latex';

% make plot pretty
xlabel('x')
ylabel('u')
latexify(19,15,20)

% mask colorbar to match transparency
annotation('rectangle',...
    bg_bar.Position,...
    'FaceAlpha',bgAlpha,...
    'EdgeColor',[1 1 1],...
    'FaceColor',[1 1 1]);

% animate and make gif
axis equal tight manual % this ensures that getframe() returns a consistent size
filename = 'test.gif';
lim_x = [a;b];
lim_y = [min(waveData(:));max(waveData(:))];
lim_z = [0;1];
update_view(lim_x,lim_y,lim_z);
playtime = T;
pausetime = playtime/size(waveData,2);

% constrain GIF to 60FPS
frames = playtime * 60;
shutter = ceil(size(waveData,2)/frames);

if play
    for i = 1:size(waveData,2)

        % update plot
        if mod(i-1,shutter)==0
            set(wave,'XData',xx,'YData',waveData(:,i))
            pause(pausetime) % pause for proper MATLAB display speed
        end
        
        % save GIF
        if saveGIF && mod(i-1,shutter)==0

            % capture the plot as an image 
            frame = getframe;
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);

            % write to the GIF file 
            if i == 1
              imwrite(imind,cm,filename,...
                        'gif','Loopcount',inf,'DelayTime',1/60);
            else
              imwrite(imind,cm,filename,...
                        'gif','WriteMode','append','DelayTime',1/60);
            end

        end

    end
end
    
hold off

%% save screenshots

sc = figure(3);

% take a fixed number of screenshots in equispaced time
num_pics = 4;
shutter = floor(size(waveData,2)/(num_pics-1));
idx = 1;

% extract the highest resolution data
waveData = uvect{end};
% create x-axis, mapped to highest resolution data
xx = linspace(a,b,size(waveData,1));

% make background gradient
bgAlpha = 0.5;
yy = linspace(min(waveData(:)),max(waveData(:)),100);
[X,Y] = meshgrid(xx,yy);
Z = c(X);

% make figure pretty
latexify(19,16)

for i = 1:num_pics

    subplot(4,1,i)
    
    % plot function
    pp = plot(xx,waveData(:,idx),'k','LineWidth',0.75);
    idx = idx + shutter;
    hold on
    
    % plot colorbar
    bg_bar = colorbar;
    bg_bar.TickLabelInterpreter = 'latex';
    bg_bar.Label.String = 'Speed';
    bg_bar.Label.Interpreter = 'latex';
    
    % mask colorbar to match transparency
    annotation('rectangle',...
        bg_bar.Position,...
        'FaceAlpha',bgAlpha,...
        'EdgeColor',[1 1 1],...
        'FaceColor',[1 1 1]);
    
    % plot background
    cmap = spring(100);
    colormap(cmap(30:60,:))
    bg = image(xx,yy,Z,'CDataMapping','scaled');
    bg.AlphaData = bgAlpha;
    
    % define axis limits
    lim_x = [a;b];
    lim_y = [min(waveData(:));max(waveData(:))];
    lim_z = [0;1];
    update_view(lim_x,lim_y,lim_z);
    
    grid on
    hold off
    
    % move wave to top of display stack
    uistack(pp,'top')
    
end

hold on
handle = axes(sc,'visible','off');
handle.XLabel.Visible = 'on';
handle.YLabel.Visible = 'on';
xlabel('x')
ylabel('u')
hold off


%% supporting functions

function update_view(lim_x,lim_y,lim_z)

    xlim(lim_x)
    ylim(lim_y)
    zlim(lim_z)

end