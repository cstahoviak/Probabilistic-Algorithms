function [ init_distr_hat, trans_prob_hat, obs_prob_hat ] = ...
    baumwelch_Mstep( p, eln_gamma, gamma, eln_xi, xi, y_obs )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if (size(gamma,2) ~= size(xi,3)+1) && (size(gamma,3) ~= size(xi,4))
    error('gamma, xi size mismatch')
    
else
    % Implement equations from Berkeley paper
%     init_distr_hat = updateInitDistr( gamma );
%     trans_prob_hat = updateTransitionProb( gamma, xi);
%     obs_prob_hat   = updateEmissionProb( p, gamma, y_obs );
    
    % Implement extended-log versions of the parameter updates
    % from the Berkeley paper:
    % WORKING!! - Now need to update functions above to achieve similar
    % results!
    init_distr_hat = updateInitDistr_eln( eln_gamma );
    trans_prob_hat = updateTransitionProb_eln( eln_gamma, eln_xi);
    obs_prob_hat   = updateEmissionProb_eln( p, eln_gamma, y_obs );
    
end

