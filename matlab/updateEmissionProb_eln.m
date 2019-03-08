function [ obs_prob_hat ] = updateEmissionProb_eln( p, eln_gamma, eln_xi, y_obs )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

n = size(eln_gamma,1);
T = size(eln_gamma,2);

% initialization
obs_prob_hat = zeros(p,n);

for i=1:p
    for j=1:n
        numerator = NaN;
        denominator = NaN;
        
        for k=1:T
            if y_obs(k) == i
                numerator = elnsum(numerator, eln_gamma(j,k));
            end
            denominator = elnsum(denominator, eln_gamma(j,k));
        end
        obs_prob_hat(i,j) = eexp(elnprod(numerator, -denominator));
    end
end

