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
title({'Maximum Likelihood Gaussian Parameter Estimation', ...
    'ThreeClass\_log 2D Data'},'Interpreter','latex');
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
hold on;

for i=1:3
    figure(i+1);
    ax(i+1) = gca;
    title({'Expectation Maximization for GMMs', ...
        strcat('Data Label','$\;$',num2str(i))},'Interpreter','latex');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    hold on;
    pause on;
end

figure(5);
ax(5) = gca;
title({'Expectation Maximization for GMMs', ...
    'ThreeClass\_log 2D Data'},'Interpreter','latex');
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
hold on;
pause on;

%% Maximum Likelihood Estimation - Multivariate Gaussian

% For Model A, the labels {pi_j = 1,2,3} associated with each data point
% are known. In this case, it is straightforward to compute the guassian
% parameters theta = {mu_j, sigma_j} =  {mu1, sigma1, mu2, sigma2, mu3,
% sigma3} via Maximum Likelihood Estimation (MLE).

% get one-pass ML estimates of mean and variance of each distribution

mu = zeros(2,max(y_obs(:,1)));
sigma = zeros(2,2,max(y_obs(:,1)));

idx_one = (y_obs(:,1) == 1) ;
idx_two = (y_obs(:,1) == 2);
idx_three = (y_obs(:,1) == 3);

for i=1:length(ax)
    if (i == 1) || (i == 5)
        h1 = scatter(ax(i),y_obs(idx_one,2),y_obs(idx_one,3), ...
            5, colors(y_obs(idx_one,1),:),'filled');

        h2 = scatter(ax(i),y_obs(idx_two,2),y_obs(idx_two,3), ...
            5, colors(y_obs(idx_two,1),:),'filled');

        h3 = scatter(ax(i),y_obs(idx_three,2),y_obs(idx_three,3), ...
            5, colors(y_obs(idx_three,1),:),'filled');

        hdl = legend(ax(i),[h1(1), h2(1), h3(1)],'Label 1','Label 2','Label 3');
        set(hdl,'Interpreter','latex','Location','Northwest','Autoupdate','off');
    end
end

[mu(:,1),sigma(:,:,1)] = gaussian_mle( y_obs(idx_one,2:3)' );
plot_gaussian_ellipsoid(mu(:,1), sigma(:,:,1), colors(1,:), 1, [], ax(1));

[mu(:,2),sigma(:,:,2)] = gaussian_mle( y_obs(idx_two,2:3)' );
plot_gaussian_ellipsoid(mu(:,2), sigma(:,:,2), colors(2,:), 1, [], ax(1));

[mu(:,3),sigma(:,:,3)] = gaussian_mle( y_obs(idx_three,2:3)' );
plot_gaussian_ellipsoid(mu(:,3), sigma(:,:,3), colors(3,:), 1, [], ax(1));

%% Expectation Maximization - Gaussian Mixture Models (GMMs)

% For Model B, the data labels {pi_j = 1,2,3} are known, but the variable
% Z_i, the sub-distribution labels {m = 1,2} are unknown (latent
% variables), and we can use Expectation Maximization (EM) to learn the
% parameters of the GMM - {w1_j, w2_j, mu1_j, sigma1_j, mu2_l, sigma2_j}.

% Since the data labels {pi_j = 1,2,3} are known, we can use EM to
% determine the paramter set {w1_j, w2_j, mu1_j, sigma1_j, mu2_l, sigma2_j}
% associated with each label j. This process for a given label j is
% independent of every other label j.

n = 2;      % dimensionality of data
M = 2;      % number of sub-distributions within label j      

% number of times EM algorithm will be done
N_iter = 100; 

for j=1:3
    
    % get data with label j
    idx = (y_obs(:,1) == j);
    X = y_obs(idx,2:3)';
    N_points = size(X,2);
    
    % init parameters
    mu      = zeros(n,N_iter,M);
    sigma   = zeros(n,n,N_iter,M);
    weights = zeros(M,N_iter);
    r       = zeros(M,N_points,N_iter);
    
    % get intial parameter estimates
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
        % clear axes and replot data
        if s < N_iter
            cla(ax(j+1));
        end
        scatter(ax(j+1),y_obs(idx,2),y_obs(idx,3), 5, ...
                colors(j,:),'filled');

        % Gaussian EM - E-step
        [ r(:,:,s) ] = gaussian_em_estep( X, squeeze(mu(:,s,:)), ...
            squeeze(sigma(:,:,s,:)), weights(:,s) );

        % Guassian EM - M-step
        [ mu(:,s+1,:), sigma(:,:,s+1,:), weights(:,s+1) ] = ...
            gaussian_em_mstep( X, r(:,:,s), squeeze(mu(:,s,:)) );

        for m=1:M
            % plot guassian ellipsoid m within label j ot sub-figure
            plot_gaussian_ellipsoid(mu(:,s+1,m), sigma(:,:,s+1,m), ...
                colors(j,:), 1, [], ax(j+1));
        end
        pause(.1);
        
    end
    
    % plot m sub-distributions within label j to main figure
    for m=1:M
        % plot guassian ellipsoid m
        plot_gaussian_ellipsoid(mu(:,s+1,m), sigma(:,:,s+1,m), ...
            colors(j,:), 1, [], ax(5));
    end
    
    fprintf('data label %d:\n', j);
    disp([mu(:,end,1), mu(:,end,2)])
    disp(sigma(:,:,end,1))
    disp(sigma(:,:,end,2))
    disp(weights(:,end)); 
    
end