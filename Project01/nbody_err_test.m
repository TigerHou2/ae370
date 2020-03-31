%% nbody_err_test.m
%
% Author: Tiger Hou
% 
% Description:
%   Runs n-body simulations to determine the error reduction rate for
%   local interpolation numerical methods such as RK4.
%
%% Setup

% clear console and workspace
close all hidden
clear;clc

% load n-body case
case_select = 6;
sub_case = 11;
nbody_test_cases;

% order reduction for gravity model
rv0 = [r0;v0];

% define various timestep resolutions
dtvect = [dt*2^7, dt*2^6, dt*2^5, dt*2^4, dt*2^3, dt*2^2, dt*2, dt];


%% Simulation

rv_end = zeros(6,size(rv0,2),length(dtvect));

for j = 1:length(dtvect)
    
    rv_end(:,:,j) = rk4(@gravity,rv0,T,dtvect(j),m,G);
    
end


%% Error Plots

errs = zeros(size(rv0,2),length(dtvect)-1);

for j = 1:length(dtvect)-1
    
    for k = 1:size(rv0,2)
        
        errs(k,j) = norm(rv_end(1:3,k,j)-rv_end(1:3,k,end)) / ...
                    norm(rv_end(1:3,k,end));
    
    end
    
end

for k = 1:size(rv0,2)
    
    figure(k)
    loglog(dtvect(1:end-1),errs(k,:))
    
end

