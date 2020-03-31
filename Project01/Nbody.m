%% Simulation

close all hidden
clear;clc

% display settings
drawsize = 20000;
ref_frame = 1;
rot_frame = [1,2];
maxNumPoints = 1e6; % max: 1e6

% load n-body case
case_select = 3;
sub_case = 11;
nbody_test_cases;

num_steps = floor(T/dt);

drawsize = min(drawsize,num_steps);

r_bary = zeros(3,1);
for j = 1:size(r0,2)
    r_bary = r_bary + m(j) * r0(:,j);
end
r_bary = r_bary / sum(m);
v_bary = zeros(3,1);
for j = 1:size(r0,2)
    v_bary = v_bary + m(j) * v0(:,j);
end
v_bary = v_bary / sum(m);
predicted_barycenter = r_bary + (T-tk)*v_bary;

% initialize iterates
r_k   = r0;
v_k   = v0;
r_km1 = r0;
v_km1 = v0;
r_kp1 = r0;
v_kp1 = v0;

rv_k = [r0;v0];
rv_kp1 = rv_k;

% advance to final time T
r_hist = zeros(3,size(r0,2),num_steps);

tic

for jj = 1 : num_steps
    
    % ============= HEU =============
%     if jj < 2
%         % propagate using Forward Euler on first step
%         v_kp1 = v_k + dt * f(r_k,m,G);
%         r_kp1 = r_k + dt * v_kp1;
%     else
%         % propagate using Heun's method
%         v_kp1 = v_k + dt * f(r_k,m,G);
%         r_kp1 = r_k + dt/2 * (v_k + v_km1)...
%                     + dt^2/2 * (f(r_km1,m,G) + f(r_k,m,G));
%     end
    
    % ============= RK4 Accel + Forward Euler Velocity =============
%     v_kp1 = v_k + dt * f(r_k,m,G);
%     k1 = dt * v_k + dt^2 * f( r_k,m,G);
%     k2 = dt * v_k + dt^2 * f((r_k + k1/2),m,G);
%     k3 = dt * v_k + dt^2 * f((r_k + k2/2),m,G);
%     k4 = dt * v_k + dt^2 * f((r_k + k3)  ,m,G);
%     r_kp1 = r_k + 1/6 * (k1 + 2*k2 + 2*k3 + k4);

    % update iterates
%     r_km1 = r_k;
%     v_km1 = v_k;
%     r_k = r_kp1;
%     v_k = v_kp1;
%     r_hist(:,:,jj) = r_kp1;
    
    % ============= RK4 =============
    k1 = dt * ff( rv_k,m,G);
    k2 = dt * ff((rv_k + k1/2),m,G);
    k3 = dt * ff((rv_k + k2/2),m,G);
    k4 = dt * ff((rv_k + k3)  ,m,G);
    rv_kp1 = rv_k + 1/6 * (k1 + 2*k2 + 2*k3 + k4);

    % update iterates
    rv_k = rv_kp1;
    r_hist(:,:,jj) = rv_kp1(1:3,:);
    
end

disp('Simulation complete!')

toc

filename = strcat('Data\Case_',num2str(case_select));
save(filename, ...
     'r_hist', 'tk', 'T', 'dt', ...
     'drawsize', 'ref_frame', 'maxNumPoints', 'num_steps', 'lgd', ...
     'case_select', 'r0', 'v0', 'G', 'm', ...
     'predicted_barycenter');

%% Animation

% create tracer line and particles
     tracer(1:size(r0,2)) = animatedline;
tracer_part(1:size(r0,2)) = animatedline;
tracer_bary(1:size(r0,2)) = animatedline;
 tracer_rot(1:size(r0,2)) = animatedline;
       anim(1:size(r0,2)) = animatedline;
  anim_part(1:size(r0,2)) = animatedline;
  anim_bary(1:size(r0,2)) = animatedline;
   anim_rot(1:size(r0,2)) = animatedline;
cmap = hsv(size(r0,2));

% ==================== Animate in Global Frame ====================

% set up animated display
figure(1)
title('Global Frame')

for i = 1:size(r0,2)
    anim(i) = animatedline(gca, 'Color', cmap(i,:), 'LineWidth', 1.5, ...
        'MaximumNumPoints', maxNumPoints);
    tracer(i) = animatedline(gca, 'MaximumNumPoints', 1, 'Marker', 'o', ...
        'MarkerEdgeColor', cmap(i,:), 'MarkerFaceColor', cmap(i,:), ...
        'Color', cmap(i,:));
end

pbaspect([1,1,1])
axis equal
[lim_x,lim_y,lim_z] = init_view(r_hist);

for j = 1:drawsize:num_steps-drawsize+1
    for i = 1:size(r0,2)
        x = r_hist(1,i,j:j+drawsize-1);
        y = r_hist(2,i,j:j+drawsize-1);
        z = r_hist(3,i,j:j+drawsize-1);
        addpoints(anim(i),x(:),y(:),z(:))
        addpoints(tracer(i),x(end),y(end),z(end))
        [lim_x,lim_y,lim_z] = set_view_bounds(lim_x,lim_y,lim_z,x,y,z);
    end
    update_view(lim_x,lim_y,lim_z);
    legend(anim,lgd)
    drawnow
end


% ==================== Animate in Particle Frame ====================

% change reference frame and animate
% ref_frame = 2; % this is now define at the top of the script
figure(2);
title([lgd(ref_frame) ' Frame'])

% set up animated display
for i = 1:size(r0,2)
    anim_part(i) = animatedline(gca, 'Color', cmap(i,:), 'LineWidth', 1.5, ...
        'MaximumNumPoints', maxNumPoints);
    tracer_part(i) = animatedline(gca, 'MaximumNumPoints', 1, 'Marker', 'o', ...
        'MarkerEdgeColor', cmap(i,:), 'MarkerFaceColor', cmap(i,:), ...
        'Color', cmap(i,:));
end

r_hist_frame = r_hist;
for i = 1:num_steps
    r_hist_frame(:,:,i) = r_hist_frame(:,:,i) - r_hist_frame(:,ref_frame,i);
end

pbaspect([1,1,1])
axis equal
[lim_x,lim_y,lim_z] = init_view(r_hist_frame);

for j = 1:drawsize:num_steps-drawsize+1
    for i = 1:size(r0,2)
        x = r_hist_frame(1,i,j:j+drawsize-1);
        y = r_hist_frame(2,i,j:j+drawsize-1);
        z = r_hist_frame(3,i,j:j+drawsize-1);
        addpoints(anim_part(i),x(:),y(:),z(:))
        addpoints(tracer_part(i),x(end),y(end),z(end))
        [lim_x,lim_y,lim_z] = set_view_bounds(lim_x,lim_y,lim_z,x,y,z);
    end
    update_view(lim_x,lim_y,lim_z);
    legend(anim_part,lgd)
    drawnow
end


% ==================== Animate in Barycenter Frame ====================

% change reference frame and animate
figure(3);
title('Barycenter Frame')

% set up animated display
for i = 1:size(r0,2)
    anim_bary(i) = animatedline(gca, 'Color', cmap(i,:), 'LineWidth', 1.5, ...
        'MaximumNumPoints', maxNumPoints);
    tracer_bary(i) = animatedline(gca, 'MaximumNumPoints', 1, 'Marker', 'o', ...
        'MarkerEdgeColor', cmap(i,:), 'MarkerFaceColor', cmap(i,:), ...
        'Color', cmap(i,:));
end

r_hist_bary = r_hist;
for i = 1:num_steps
    barycenter = zeros(3,1);
    for j = 1:size(r0,2)
        barycenter = barycenter + m(j) * r_hist_bary(:,j,i);
    end
    barycenter = barycenter / sum(m);
    r_hist_bary(:,:,i) = r_hist_bary(:,:,i) - barycenter;
end

pbaspect([1,1,1])
axis equal
[lim_x,lim_y,lim_z] = init_view(r_hist_bary);

for j = 1:drawsize:num_steps-drawsize+1
    for i = 1:size(r0,2)
        x = r_hist_bary(1,i,j:j+drawsize-1);
        y = r_hist_bary(2,i,j:j+drawsize-1);
        z = r_hist_bary(3,i,j:j+drawsize-1);
        addpoints(anim_bary(i),x(:),y(:),z(:))
        addpoints(tracer_bary(i),x(end),y(end),z(end))
        [lim_x,lim_y,lim_z] = set_view_bounds(lim_x,lim_y,lim_z,x,y,z);
    end
    update_view(lim_x,lim_y,lim_z);
    legend(anim_bary,lgd)
    drawnow
end

disp(['Barycenter Error:' num2str(norm(barycenter-predicted_barycenter))])


% ==================== Animate in Rotating Frame ====================

% change reference frame and animate
figure(4);
title(['Rotating Frame of ' lgd{rot_frame(1)} ' and ' lgd{rot_frame(2)}])

% set up animated display
for i = 1:size(r0,2)
    anim_rot(i) = animatedline(gca, 'Color', cmap(i,:), 'LineWidth', 1.5, ...
        'MaximumNumPoints', maxNumPoints);
    tracer_rot(i) = animatedline(gca, 'MaximumNumPoints', 1, 'Marker', 'o', ...
        'MarkerEdgeColor', cmap(i,:), 'MarkerFaceColor', cmap(i,:), ...
        'Color', cmap(i,:));
end

r_hist_rot = r_hist;
for i = 1:num_steps
    translation = -r_hist_rot(:,rot_frame(1),i);
    cur_vect = r_hist_rot(:,rot_frame(2),i) - r_hist_rot(:,rot_frame(1),i);
    tgt_vect = [1,0,0];
    rotation = vrrotvec2mat(vrrotvec(cur_vect(:),tgt_vect(:)));
    for j = 1:size(r0,2)
        r_hist_rot(:,j,i) = r_hist_rot(:,j,i) + translation;
        r_hist_rot(:,j,i) = rotation * reshape(r_hist_rot(:,j,i),[3,1]);
    end
end

pbaspect([1,1,1])
axis equal
[lim_x,lim_y,lim_z] = init_view(r_hist_rot);

for j = 1:drawsize:num_steps-drawsize+1
    for i = 1:size(r0,2)
        x = r_hist_rot(1,i,j:j+drawsize-1);
        y = r_hist_rot(2,i,j:j+drawsize-1);
        z = r_hist_rot(3,i,j:j+drawsize-1);
        addpoints(anim_rot(i),x(:),y(:),z(:))
        addpoints(tracer_rot(i),x(end),y(end),z(end))
        [lim_x,lim_y,lim_z] = set_view_bounds(lim_x,lim_y,lim_z,x,y,z);
    end
    update_view(lim_x,lim_y,lim_z);
    legend(anim_rot,lgd)
    drawnow
end



%% Function Definitions

function a = f(r,m,G)

N = size(r,2);
a = zeros(size(r));

for i = 1:N
    for j = 1:N
        if i ~= j
        a(:,i) = a(:,i) + G*m(j)*(r(:,j)-r(:,i)) / norm(r(:,j)-r(:,i))^3;
        end
    end
end

end

function rv_dot = ff(rv,m,G)

N = size(rv,2);
rv_dot = zeros(size(rv));

for i = 1:N
    % velocity
    rv_dot(1:3,i) = rv(4:6,i);
    % acceleration
    for j = 1:N
        if i ~= j
        rv_dot(4:6,i) = rv_dot(4:6,i) + ...
                        G*m(j)*(rv(1:3,j)-rv(1:3,i)) / norm(rv(1:3,j)-rv(1:3,i))^3;
        end
    end
end

end

function [lim_x,lim_y,lim_z] = init_view(trajectory)

    lim_x = [min(trajectory(1,:,1)), max(trajectory(1,:,1))];
    if (lim_x(1) == lim_x(2)) 
        lim_x = [-1,1]; end
    lim_y = [min(trajectory(2,:,1)), max(trajectory(2,:,1))];
    if (lim_y(1) == lim_y(2)) 
        lim_y = [-1,1]; end
    lim_z = [min(trajectory(3,:,1)), max(trajectory(3,:,1))];
    if (lim_z(1) == lim_z(2)) 
        lim_z = [-1,1]; end
    
    update_view(lim_x,lim_y,lim_z);

end

function [lim_x,lim_y,lim_z] = set_view_bounds(lim_x,lim_y,lim_z,x,y,z)

    lim_x = [min([lim_x(1);x(:)]), max([lim_x(2);x(:)])];
    lim_y = [min([lim_y(1);y(:)]), max([lim_y(2);y(:)])];
    lim_z = [min([lim_z(1);z(:)]), max([lim_z(2);z(:)])];

end

function update_view(lim_x,lim_y,lim_z)

    xlim(lim_x)
    ylim(lim_y)
    zlim(lim_z)

end