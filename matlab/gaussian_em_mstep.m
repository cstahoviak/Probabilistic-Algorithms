function [ mu, sigma, weights ] = ...
    gaussian_em_mstep( X, r, mu )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

n = size(mu,1);         % dimensionality of data
M = size(mu,2);         % number of sub-distributions within label j
Npoints = size(X,2);    % number of data points in label j

% fprintf('mstep: [M, Npoints] = [%d, %d]\n', M, Npoints);

% initialize
mu_new = zeros(n,M);
sigma_new = zeros(n,n,M);
weights_new = zeros(M,1);

% Gaussian EM - M-step
for m=1:M
    for i=1:Npoints
        % update mean of distribution m
        mu_new(:,m) = mu_new(:,m) + r(m,i)*X(:,i);

        % update covariance matrix of distribution m
        sigma_new(:,:,m) = sigma_new(:,:,m) + r(m,i) * ...
            (X(:,i) - mu(:,m))*(X(:,i) - mu(:,m))';

        % update weight of distribution m
        weights_new(m,1) = weights_new(m,1) + r(m,i);
    end

    % sanity check - do these match?
%     fprintf('Nm(%d) = %f\n', m, weights_new(m,1))
%     fprintf('Nm(%d) = sum(r(m,:)) =  %f\n\n', m, sum(r(m,:)))

    % current number of points assigned to cluster m
    Nm = sum(r(m,:));

    % normlize by total number of data points
    mu_new(:,m)      = mu_new(:,m)/Nm;
    sigma_new(:,:,m) = sigma_new(:,:,m)/Nm;
    weights_new(m,1) = weights_new(m,1)/Npoints;
    
end
    
    % resize parameters
    mu = reshape(mu_new,n,1,M);
    sigma = reshape(sigma_new,n,n,1,M);
    weights = weights_new;
    
end

