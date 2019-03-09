function [ eln_alpha, eln_beta, eln_gamma, gamma, eln_xi, xi ] = ...
    baumwelch_Estep( init_distr, trans_prob, obs_prob, y_obs )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Mann log-weighted (numerically-stable) forward-backward alg.
eln_alpha = forward_eln( init_distr, trans_prob, obs_prob, y_obs );
eln_beta = backward_eln( trans_prob, obs_prob, y_obs );

% get the log-posterior and posterior distributions (exact inference)
% [ eln_gamma, gamma ] = posterior_elnfb( eln_alpha(:,2:end), eln_beta );
[ eln_gamma, gamma ] = posterior_elnfb( eln_alpha, eln_beta );

% get log of probability xi(i,j,k)
% [ eln_xi, xi ] = elnxi( eln_alpha(:,2:end), eln_beta(:,2:end), ...
%     trans_prob, obs_prob, y_obs );

[ eln_xi, xi ] = elnxi( eln_alpha, eln_beta, ...
    trans_prob, obs_prob, y_obs );

end

