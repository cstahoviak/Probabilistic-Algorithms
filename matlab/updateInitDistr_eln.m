function [ init_distr_hat ] = updateInitDistr_eln( eln_gamma )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

n = size(eln_gamma,1);
D = size(eln_gamma,3);

init_distr_hat = zeros(n,1);

% get initial distribution estimate
for d=1:D 
    for i=1:n
        init_distr_hat(i,1) = init_distr_hat(i,1) + eexp(eln_gamma(i,1,d));
    end
end
    
init_distr_hat = init_distr_hat/D;

% normalize
% init_distr_hat = init_distr_hat./sum(init_distr_hat);

end