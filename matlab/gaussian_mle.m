function [ mu, sigma ] = gaussian_mle( X )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n = size(X,1);  % number of variables in multivariate Gaussian
M = size(X,2);  % number of observations of each n-dim realization of X

mu = sum(X,2)/M;

sigma = zeros(n);
for i=1:size(X,2)
    sigma = sigma + (X(:,i) - mu)*(X(:,i) - mu)';
end

sigma = sigma/M;

