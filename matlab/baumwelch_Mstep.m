function [ init_distr_hat, trans_prob_hat, obs_prob_hat ] = ...
    baumwelch_Mstep( p, gamma, xi, y_obs )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if (size(gamma,2) ~= size(xi,3)+1) && (size(gamma,3) ~= size(xi,4))
    error('gamma, xi size mismatch')
    
else
    
    % Implement equations from Berkeley paper
    init_distr_hat = updateInitDistr( gamma );
    trans_prob_hat = updateTransitionProb( gamma(:,2:end), xi);
    obs_prob_hat   = upsateEmissionProb( p, gamma(:,2:end), y_obs );

end

