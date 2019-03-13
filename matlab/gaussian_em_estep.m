function [ r ] = gaussian_em_estep( X, mu, sigma, weights )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

M = size(mu,2);         % number of sub-distributions within label j
Npoints = size(X,2);    % number of data points in label j

% fprintf('estep: [M, Npoints] = [%d, %d]\n', M, Npoints);  

% init responsibility
r = zeros(M,Npoints);

% Gaussian EM - E-step
for m = 1:M
    for i = 1:Npoints
        numerator = weights(m,1) * ...
            mvnpdf(X(:,i),mu(:,m),sigma(:,:,m));

        denominator = 0;
        for j=1:M
            denominator = denominator + weights(j,1) * ...
                mvnpdf(X(:,i),mu(:,j),sigma(:,:,j));
        end

        % compute "responsibility" that data point i belongs to cluster m
        r(m,i) = numerator/denominator;
    end
end

end

