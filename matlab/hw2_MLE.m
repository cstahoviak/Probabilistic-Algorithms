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

%% Plot Data

colors = [0,      0.4470, 0.7410;
          0.9290, 0.6940, 0.1250;
          0.6350, 0.0780, 0.1840];

figure(1);
ax1 = gca;
title({'Maximum Likelihood Gaussian Parameter Estimation', ...
    'ThreeClass\_log 2D Data'},'Interpreter','latex');
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
hold on;

figure(2);
ax2 = gca;
title({'Expectation Maximization for GMMs', ...
    'ThreeClass\_log 2D Data'},'Interpreter','latex');
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
hold on;

for i=1:size(y_obs,1)
    % maximum likelihood plot
    scatter(ax1,y_obs(i,2),y_obs(i,3), 5, ...
        colors(y_obs(i,1),:),'filled');
    
    % expectation maximization plot
    scatter(ax2,y_obs(i,2),y_obs(i,3), 5, ...
        colors(y_obs(i,1),:),'filled');
end

%% Maximum Likelihood Estimation - Multivariate Gaussian

% For Model A, the labels {pi_j = 1,2,3} associated with each data point
% are known. In this case, it is straightforward to compute the guassian
% parameters theta = {mu_j, sigma_j} =  {mu1, sigma1, mu2, sigma2, mu3,
% sigma3} via Maximum Likelihood Estimation (MLE).

% get one-pass ML estimates of mean and variance of each distribution

idx_one = (y_obs(:,1) == 1);
idx_two = (y_obs(:,1) == 2);
idx_three = (y_obs(:,1) == 3);

[mu1,sigma1] = gaussian_mle( y_obs(idx_one,2:3)' );
plot_gaussian_ellipsoid(mu1, sigma1, colors(1,:), 1, [], ax1);

[mu2,sigma2] = gaussian_mle( y_obs(idx_two,2:3)' );
plot_gaussian_ellipsoid(mu2, sigma2, colors(2,:), 1, [], ax1);

[mu3,sigma3] = gaussian_mle( y_obs(idx_three,2:3)' );
plot_gaussian_ellipsoid(mu3, sigma3, colors(3,:), 1, [], ax1);

%% Expectation Maximization - Gaussian Mixture Models (GMMs)

% For Model B, the data labels {pi_j = 1,2,3} are known, but the variable
% Z_i, the sub-distribution labels {m = 1,2} are unknown (latent
% variables), and we can use Expectation Maximization (EM) to learn the
% parameters of the GMM - {w1_j, w2_j, mu1_j, sigma1_j, mu2_l, sigma2_j}.

% Since the data labels {pi_j = 1,2,3} are known, we can use EM to
% determine the paramter set {w1_j, w2_j, mu1_j, sigma1_j, mu2_l, sigma2_j}
% associated with each label j. This process for a given label j is
% independent of every other label j.









