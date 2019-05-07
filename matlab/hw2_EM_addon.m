%% Header

% Filename:     hw2_MLE.m
% Author:       Carl Stahoviak
% Date Created: 03/02/2019

clc;
clear;
close ALL;

%% Load Data

% load MLE Data Log
load('ThreeClass_log.mat')
y_obs = ThreeClass_log;

%% Initialize Data Plots

colors = [0,      0.4470, 0.7410;
          0.9290, 0.6940, 0.1250;
          0.6350, 0.0780, 0.1840;
          0.4940, 0.1840, 0.5560;
          0.4660, 0.6740, 0.1880];

figure(1);
ax(1) = gca;
title({'Expectation Maximization for GMMs', ...
    'Unknown Data Labels - Colored by Truth Label'},'Interpreter','latex');
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
hold on;
pause on;

figure(2);
ax(2) = gca;
title({'Expectation Maximization for GMMs', ...
    'Unknown Data Labels - Colored by Discovered Label'},'Interpreter','latex');
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
hold on;
pause on;

%% Expectation Maximization - Gausssian Mixture Models (GMMs) - Part 2

% In this section, we will assume that the data labels j={1,2,3} are
% unknown, and thus the Incomplete Log-Likelihood will be solved via the
% Expected Complete Log-Likelihood using EM.

n = 2;      % dimensionality of data
M = 3;      % number of sub-distributions within label j      

% number of times EM algorithm will be done
N_iter = 100; 

X = y_obs(:,2:3)';
N_points = size(X,2);

% init parameters
mu      = zeros(n,N_iter,M);
sigma   = zeros(n,n,N_iter,M);
weights = zeros(M,N_iter);
r       = zeros(M,N_points,N_iter);
labels  = zeros(N_points,N_iter);

% get intial parameter estimates
% NOTE: interesting initialization given by:
% mu = [3.9158   -5.4860   -4.3932;
%      -1.8518   15.2275    6.5673];

lower = [min(X(1,:)), min(X(2,:))];   % min [x,y]
upper = [max(X(1,:)), max(X(2,:))];   % max [x,y]
for m=1:M
    for i=1:n
        mu(i,1,m) = (upper(i)-lower(i)).*rand(1) + lower(i);
    end
    sigma(:,:,1,m) = [5, 0.5; 0.5 5];
end
weights(:,1) = ones(M,1)/M;

for s=1:N_iter
    % Gaussian EM - E-step
    [ r(:,:,s) ] = gaussian_em_estep( X, squeeze(mu(:,s,:)), ...
        squeeze(sigma(:,:,s,:)), weights(:,s) );
    
    % generate updated labels for each data point i based on r(m,i)
    [~,labels(:,s)] = max(r(:,:,s));

    % Guassian EM - M-step
    [ mu(:,s+1,:), sigma(:,:,s+1,:), weights(:,s+1) ] = ...
        gaussian_em_mstep( X, r(:,:,s), squeeze(mu(:,s,:)) );
    
    % clear axes and replot data
    cla(ax(1));
    cla(ax(2));
    % plot data with truth labels
    scatter(ax(1), y_obs(:,2), y_obs(:,3), 5, ...
        colors(y_obs(:,1),:),'filled');
    % plot data with discovered labels
    scatter(ax(2), y_obs(:,2), y_obs(:,3), 5, ...
        colors(labels(:,s),:),'filled');
    
%     scatter(ax(i),y_obs(idx_one,2),y_obs(idx_one,3), ...
%             5, colors(y_obs(idx_one,1),:),'filled');

    for m=1:M
        % plot guassian ellipsoid m
        plot_gaussian_ellipsoid(mu(:,s+1,m), sigma(:,:,s+1,m), ...
            'k', 1, [], ax(1));
        plot_gaussian_ellipsoid(mu(:,s+1,m), sigma(:,:,s+1,m), ...
            'k', 1, [], ax(2));
    end
    pause(.1);
end

