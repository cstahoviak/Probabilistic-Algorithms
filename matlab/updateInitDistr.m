function [ init_distr_hat ] = updateInitDistr( gamma )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

n = size(gamma,1);
D = size(gamma,3);

init_distr_hat = zeros(n,1);

% get initial distribution estimate
for i=1:n
    init_distr_hat(i,1) = sum(gamma(i,1,:))/D;
end

% normalize
init_distr_hat = init_distr_hat./sum(init_distr_hat);

end

