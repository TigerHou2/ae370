%% Problem 2 Part a
close all
clear;clc

%--problem parameters

    b = 2; %birth rate of prey
    p = 1; %effect of predation on prey
    d = 3; %death rate of predators
    g = 1; %growth of predators due to eating prey
    
    t0 = 0; %starting time
    T = 50; %final time
    
    u0 = [1;1]; %IC
    
    f = @(u) [(b-p*u(2))*u(1); (g*u(1)-d)*u(2)]; %RHS
    
%--

%part a)
    %--simulation params
        dtvect = [5e-2, 2.5e-2, 1e-2, 5e-3, 2.5e-3, 1e-3, 5e-4, 2.5e-4]; 
    %--

    %initialize vector that stores approximate soln at T for various dt
    u_FWE_keep = zeros( 2,length( dtvect ) ); % forward Euler
    u_AB2_keep = zeros( 2,length( dtvect ) ); % Adams-Bashforth 2-step
    u_HEU_keep = zeros( 2,length( dtvect ) ); % Heun's
    u_RK4_keep = zeros( 2,length( dtvect ) ); % Runge-Kutta 4
    
    
    %advance in time
    for j = 1 : length(dtvect)

        %current dt
        dt = dtvect(j);
        
        tk = t0; %initialize time iterate
        
        %initialize iterates for various methods
        % FWE
        u_FWE_k = u0;
        % AB2
        u_AB2_k = u0; 
        u_AB2_km1 = u0; 
        % HEU
        u_HEU_k = u0;
        % RK4
        u_RK4_k = u0;
       
        
        %run to final time T
        for jj = 1 : T/dt

            % ============= FWE =============
            % propagate and update
            u_FWE_k = u_FWE_k + f(u_FWE_k) * dt;
            
            % ============= AB2 =============
            % propagate
            if jj < 2
                %advance with Heun's for 1st time step
                u_AB2_kp1 = u_AB2_k + ...
                            0.5 * dt * ...
                            ( f(u_AB2_k) + f(u_AB2_k + dt * f(u_AB2_k)) );
            else
                u_AB2_kp1 = u_AB2_k + dt / 2 * ...
                            ( -f(u_AB2_km1) + 3 * f(u_AB2_k) );
            end
            % update iterates
            u_AB2_km1 = u_AB2_k;
            u_AB2_k = u_AB2_kp1;
            
            % ============= HEU =============
            % propagate and update
            u_HEU_k = u_HEU_k + ...
                      0.5 * dt * ...
                      ( f(u_HEU_k) + f(u_HEU_k + dt * f(u_HEU_k)) );
            
            % ============= RK4 =============
            % propagate and update
            y1 = f(u_RK4_k);
            y2 = f(u_RK4_k + dt/2 * y1);
            y3 = f(u_RK4_k + dt/2 * y2);
            y4 = f(u_RK4_k + dt * y3);
            u_RK4_k = u_RK4_k + 1/6 * dt * (y1 + 2*y2 + 2*y3 + y4);
            
            
            % ============= update time =============

            tk = tk + dt;
         
        end
        
        % store solution at T
        % FWE
        u_FWE_keep(:,j) = u_FWE_k;
        % AB2
        u_AB2_keep(:,j) = u_AB2_kp1;
        % HEU
        u_HEU_keep(:,j) = u_HEU_k;
        % RK4
        u_RK4_keep(:,j) = u_RK4_k;
        
    end
    
    %compute difference between solution at smallest dt and the other dts
    
    %initialize vector
    u_FWE_diff = zeros( length( dtvect )-1,1 );
    u_AB2_diff = zeros( length( dtvect )-1,1 );
    u_HEU_diff = zeros( length( dtvect )-1,1 );
    u_RK4_diff = zeros( length( dtvect )-1,1 );
    
    for j = 1 : length( dtvect )-1
        
        u_FWE_diff(j) = norm(u_FWE_keep(:,j) - u_FWE_keep(:,end)) ...
                            / norm(u_FWE_keep(:,end));
        u_AB2_diff(j) = norm(u_AB2_keep(:,j) - u_AB2_keep(:,end)) ...
                            / norm(u_AB2_keep(:,end));
        u_HEU_diff(j) = norm(u_HEU_keep(:,j) - u_HEU_keep(:,end)) ...
                            / norm(u_HEU_keep(:,end));
        u_RK4_diff(j) = norm(u_RK4_keep(:,j) - u_RK4_keep(:,end)) ...
                            / norm(u_RK4_keep(:,end));
       
    end
    
    f = figure(1);
    
    subplot(2,2,1)
    loglog( dtvect(1:end-1), u_FWE_diff, '.-', 'markersize', 25, 'linewidth', 2 )
    title('FWE', 'fontsize', 16, 'interpreter' , 'latex')
    set( gca, 'Color', [1 1 1] )
    set( gca, 'fontsize', 16, 'ticklabelinterpreter', 'latex' )
    xlabel('$\Delta t$', 'fontsize', 16, 'interpreter' , 'latex')
    grid(gca,'minor')
    grid on
    
    subplot(2,2,2)
    loglog( dtvect(1:end-1), u_AB2_diff, '.-', 'markersize', 25, 'linewidth', 2 )
    title('AB2', 'fontsize', 16, 'interpreter' , 'latex')
    set( gca, 'Color', [1 1 1] )
    set( gca, 'fontsize', 16, 'ticklabelinterpreter', 'latex' )
    xlabel('$\Delta t$', 'fontsize', 16, 'interpreter' , 'latex')
    grid(gca,'minor')
    grid on
    
    subplot(2,2,3)
    loglog( dtvect(1:end-1), u_HEU_diff, '.-', 'markersize', 25, 'linewidth', 2 )
    title('HEU', 'fontsize', 16, 'interpreter' , 'latex')
    set( gca, 'Color', [1 1 1] )
    set( gca, 'fontsize', 16, 'ticklabelinterpreter', 'latex' )
    xlabel('$\Delta t$', 'fontsize', 16, 'interpreter' , 'latex')
    grid(gca,'minor')
    grid on
    
    subplot(2,2,4)
    loglog( dtvect(1:end-1), u_RK4_diff, '.-', 'markersize', 25, 'linewidth', 2 )
    title('RK4', 'fontsize', 16, 'interpreter' , 'latex')
    set( gca, 'Color', [1 1 1] )
    set( gca, 'fontsize', 16, 'ticklabelinterpreter', 'latex' )
    xlabel('$\Delta t$', 'fontsize', 16, 'interpreter' , 'latex')
    grid(gca,'minor')
    grid on
    
    
    sgtitle('$\frac{||u_{\Delta t} - u_{2.5\times10^{-4}}||}{|| u_{2.5\times10^{-4}} ||}$', 'fontsize', 28,  'interpreter' , 'latex' )
    set(f, 'PaperPositionMode', 'manual')
    set(f, 'Color', [1 1 1])
    set(f, 'PaperUnits', 'centimeters')
    set(f, 'PaperSize', [30 30])
    set(f, 'Units', 'centimeters' )
    set(f, 'Position', [0 0 30 30])
    set(f, 'PaperPosition', [0 0 30 30])

    svnm = 'error_q2';
    print( '-dpng', svnm, '-r200' )

    % store the reference values for all four methods in a .mat file
    u_FWE_ref = u_FWE_keep(:,end);
    u_AB2_ref = u_AB2_keep(:,end);
    u_HEU_ref = u_HEU_keep(:,end);
    u_RK4_ref = u_RK4_keep(:,end);
    
    save('ref_data.mat','u_FWE_ref','u_AB2_ref','u_HEU_ref','u_RK4_ref');
 
%% Problem 2 Part b, AB2
close all
clear;clc

% maunally edit this variable for convergence
num_steps = 43233;

% simulation params
b = 2; %birth rate of prey
p = 1; %effect of predation on prey
d = 3; %death rate of predators
g = 1; %growth of predators due to eating prey
u0 = [1;1]; %IC
f = @(u) [(b-p*u(2))*u(1); (g*u(1)-d)*u(2)]; %RHS

tk = 0; % starting time
T = 50; % final time
dt = (T-tk)/num_steps; % time step

% load reference value
load('ref_data.mat','u_AB2_ref');

% ====================== BEGIN TIMER ========================
tic

%initialize iterates for AB2
u_AB2_k = u0;
u_AB2_km1 = u0;

%advance to final time T
for jj = 1 : T/dt
    
    % ============= AB2 =============
    % propagate
    if jj < 2
        %advance with Heun's for 1st time step
        u_AB2_kp1 = u_AB2_k + ...
                    0.5 * dt * ...
                    ( f(u_AB2_k) + f(u_AB2_k + dt * f(u_AB2_k)) );
    else
        u_AB2_kp1 = u_AB2_k + dt / 2 * ...
                    ( -f(u_AB2_km1) + 3 * f(u_AB2_k) );
    end
    % update iterates
    u_AB2_km1 = u_AB2_k;
    u_AB2_k = u_AB2_kp1;
    tk = tk + dt;
    
end

% ====================== END TIMER ========================
t_AB2 = toc;

% manually check for convergence
u_AB2_keep = u_AB2_kp1;
u_AB2_diff = norm(u_AB2_keep - u_AB2_ref) / norm(u_AB2_ref);
disp(['Error: ' num2str(u_AB2_diff)])
if u_AB2_diff <= 6e-5
    disp('Successfully converged!');
else
    disp('Failed to converge.');
end
disp(['Time to converge: ' num2str(t_AB2) ' seconds'])


%% Problem 2 Part b, HEU
close all
clear;clc

% maunally edit this variable for convergence
num_steps = 17362;

% simulation params
b = 2; %birth rate of prey
p = 1; %effect of predation on prey
d = 3; %death rate of predators
g = 1; %growth of predators due to eating prey
u0 = [1;1]; %IC
f = @(u) [(b-p*u(2))*u(1); (g*u(1)-d)*u(2)]; %RHS

tk = 0; % starting time
T = 50; % final time
dt = (T-tk)/num_steps; % time step

% load reference value
load('ref_data.mat','u_HEU_ref');

% ====================== BEGIN TIMER ========================
tic

%initialize iterates for HEU
u_HEU_k = u0;

%advance to final time T
for jj = 1 : T/dt
    
    % ============= HEU =============
    % propagate and update
    u_HEU_k = u_HEU_k + ...
              0.5 * dt * ...
              ( f(u_HEU_k) + f(u_HEU_k + dt * f(u_HEU_k)) );
    tk = tk + dt;
    
end

% ====================== END TIMER ========================
t_HEU = toc;

% manually check for convergence
u_HEU_keep = u_HEU_k;
u_HEU_diff = norm(u_HEU_keep - u_HEU_ref) / norm(u_HEU_ref);
disp(['Error: ' num2str(u_HEU_diff)])
if u_HEU_diff <= 6e-5
    disp('Successfully converged!');
else
    disp('Failed to converge.');
end
disp(['Time to converge: ' num2str(t_HEU) ' seconds'])


%% Problem 2 Part b, RK4
close all
clear;clc

% maunally edit this variable for convergence
num_steps = 1711;

% simulation params
b = 2; %birth rate of prey
p = 1; %effect of predation on prey
d = 3; %death rate of predators
g = 1; %growth of predators due to eating prey
u0 = [1;1]; %IC
f = @(u) [(b-p*u(2))*u(1); (g*u(1)-d)*u(2)]; %RHS

tk = 0; % starting time
T = 50; % final time
dt = (T-tk)/num_steps; % time step

% load reference value
load('ref_data.mat','u_RK4_ref');

% ====================== BEGIN TIMER ========================
tic

%initialize iterates for RK4
u_RK4_k = u0;

%advance to final time T
for jj = 1 : T/dt
    
    % ============= RK4 =============
    % propagate and update
    y1 = f(u_RK4_k);
    y2 = f(u_RK4_k + dt/2 * y1);
    y3 = f(u_RK4_k + dt/2 * y2);
    y4 = f(u_RK4_k + dt * y3);
    u_RK4_k = u_RK4_k + 1/6 * dt * (y1 + 2*y2 + 2*y3 + y4);
    tk = tk + dt;
    
end

% ====================== END TIMER ========================
t_RK4 = toc;

% manually check for convergence
u_RK4_keep = u_RK4_k;
u_RK4_diff = norm(u_RK4_keep - u_RK4_ref) / norm(u_RK4_ref);
disp(['Error: ' num2str(u_RK4_diff)])
if u_RK4_diff <= 6e-5
    disp('Successfully converged!');
else
    disp('Failed to converge.');
end
disp(['Time to converge: ' num2str(t_RK4) ' seconds'])