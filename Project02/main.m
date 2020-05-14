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

saveGIF = true;
play = true;

num_pics = 4;
saveFigs = true;

hideWallEnergy = false;

if saveFigs
    fileName = ...
    input(['Please enter the prefix for your file names:' newline],'s');
end

% define left/right boundaries
a = 0;  % left boundary
b = 8;  % right boundary
ln = b-a;

% common starting point for wall and transition point
lc = 3.0;

% define wave speed for wall cases
% tc  = 1.5; % thickness
% spA = 1;  % speed of sound not in wall
% spW = 10;  % speed of sound in wall
% c = @(x) (x< lc) * spA + ...
%          (x>=lc & x < lc+tc) * spW + ...
%          (x>=lc+tc) * spA;

% define wave speed for transition cases
spL = 1.0; % speed of sound on the left
spR = 3.0; % speed of sound on the right
sharpness = 100; % larger value = sharper transition ; use [5, 15, 100]
c = @(x) ( (spL+spR)/2 - (spL-spR)/2*atan(sharpness*(x-lc))/pi*2 );

% initial conditions
w0 = @(x) zeros(size(x));
v0 = @(x) zeros(size(x));

% boundary conditions
f = @(t) (-cos(t*2*pi)+1)*(pi/2-atan(10000*(t-1)))/pi;
g = @(t) 0;

% # of n points to use
nvect = [50, 100, 200, 500, 1000];

% final time
T = 8;
% dt must match dx in order of magnitude
dtvect = 0.9 * (b-a)./nvect/max(c(linspace(a,b,10000)));

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
        tvect = 0:dt:T-dt;
        
        % second order states initialization
        um1 = w0(xj(2:end-1));
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

loglog( ln./nvect(1:end-1), err , '.', 'markersize', 32, 'Color', [0.9 0.54 0.72])
legend('$O(\Delta x^2)$', 'Error', 'Location', 'Best');
hold off

xlabel('$\Delta x$')
ylabel('Error')

setgrid(0.3,0.9)
latexify(19,19,28)

if saveFigs
    svnm = ['figures/' fileName '_error'];
    print( '-dpng', svnm, '-r200' )
end


%% plot energy

figure(2)

% extract the highest resolution data
waveData = uvect{end};
% create x-axis, mapped to highest resolution data
xx = linspace(a,b,size(waveData,1))';

if hideWallEnergy
    % cut both waveData and xx to contain only points before boundary
          xx(ceil(size(waveData,1)*(lc-a)/ln):end,:) = [];
    waveData(ceil(size(waveData,1)*(lc-a)/ln):end,:) = [];
end

dx = ( b - a ) / nvect(end);
dt = dtvect(end);
tvect = 0:dt:T-2*dt;

indices = 1:size(waveData,2)-1;
E = nan(size(indices));
for i = 1:length(indices)
    E(i) = checkEnergy(indices(i),waveData,dx,dt,c,xx);
end
plot(tvect,E,'LineWidth',3,'Color',[0.9 0.54 0.72])

xlabel('Time')
ylabel('Energy')
setgrid(0.3,0.9)
latexify(19,10,28)
expand(0,0.04)

if hideWallEnergy
    fileExt = '_energy_prewall';
else
    fileExt = '_energy';
end
if saveFigs
    svnm = ['figures/' fileName fileExt];
    print( '-dpng', svnm, '-r200' )
end


%% plot waves

if play || saveGIF
    
figure(3)

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
filename = ['figures/' fileName '.gif'];
lim_x = [a;b];
lim_y = [min(waveData(:));max(waveData(:))];
lim_z = [0;1];
update_view(lim_x,lim_y,lim_z);
playtime = T;

% constrain GIF to 60FPS
frames = playtime * 60;
shutter = ceil(size(waveData,2)/frames);

for i = 1:size(waveData,2)

    % update plot
    if mod(i-1,shutter)==0
        set(wave,'XData',xx,'YData',waveData(:,i))
        pause(1/60) % pause for proper MATLAB display speed
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
    
hold off
    
end


%% save screenshots

sc = figure(4);

% take a fixed number of screenshots in equispaced time
shutter = floor(size(waveData,2)/num_pics);
idx = shutter;

% extract the highest resolution data
waveData = uvect{end};
% create x-axis, mapped to highest resolution data
xx = linspace(a,b,size(waveData,1));

% make background gradient
bgAlpha = 0.5;
yy = linspace(min(waveData(:)),max(waveData(:)),100);
[X,Y] = meshgrid(xx,yy);
Z = c(X);

% set axis limits
lim_x = [a;b];
lim_y = [min(waveData(:));max(waveData(:))];
lim_z = [0;1];

% make figure pretty
latexify(19,12)

for i = 1:num_pics

    subplot(4,1,i)
    
    % plot function
    pp = plot(xx,waveData(:,idx),'k','LineWidth',0.75);
    idx = idx + shutter;
    hold on
    update_view(lim_x,lim_y,lim_z);
    
    % expand to window boundary
    expand(0.05,0.15)
    
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
    
    % annotate time
    posx = lim_x(1) + 0.91 * (lim_x(2) - lim_x(1));
    posy = lim_y(1) + 0.85 * (lim_y(2) - lim_y(1));
    text(posx,posy,['t = ', num2str(T/num_pics*i)])
    
end

hold on
handle = axes(sc,'visible','off');
handle.XLabel.Visible = 'on';
handle.YLabel.Visible = 'on';
xlabel('x')
ylabel('u')
expand(0.01,0.04,0.02,0.06)

latexify(19,12,11)
    
% plot colorbar
bg_bar = colorbar;
caxis([min(c(xx))-0.01,max(c(xx))+0.01])
bg_bar.TickLabelInterpreter = 'latex';
bg_bar.Label.String = 'Speed';
bg_bar.Label.Interpreter = 'latex';

% mask colorbar to match transparency
annotation('rectangle',...
    bg_bar.Position,...
    'FaceAlpha',bgAlpha,...
    'EdgeColor',[1 1 1],...
    'FaceColor',[1 1 1]);

hold off

if saveFigs
    svnm = ['figures/' fileName '_snaps'];
    print( '-dpng', svnm, '-r200' )
end


%% safety

% if we run part of the code again, don't resave figures
saveFigs = false;


%% supporting functions

function update_view(lim_x,lim_y,lim_z)

    xlim(lim_x)
    ylim(lim_y)
    zlim(lim_z)

end

function E = checkEnergy(idx,waveData,dx,dt,c,xx)

    dudx = (waveData(2:end,idx) - waveData(1:end-1,idx)) / dx;
    dudt = (waveData(:,idx+1) - waveData(:,idx)) / dt;
    
    cc = c(xx).^2;
    
    E = trapz(xx,dudt.^2) ...
      + trapz(xx(1:end-1),cc(1:end-1).*dudx.^2);

end