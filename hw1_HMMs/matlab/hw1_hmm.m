%% Header

% Author: Carl Stahoviak
% Date Created:

clc;
clear;
close ALL;

%% Load data

load('nominal_hmm_params.mat')

trans_prob = pxk_xkm1;
obs_prob = pyk_xk;

%% The Forward-Backward Algorithm (+ Ext-Log FB Alg.)

% load('nominal_hmm_short_log.mat')
load('nominal_hmm_long_log.mat')

%%% standard forward-backward algorithm
[alpha, alpha2] = forward( px0, trans_prob, obs_prob, y_obs );    % Txn
[beta, beta2] = backward( trans_prob, obs_prob, y_obs );                    % Txn
% get the posterior distribution
posterior = fb_posterior( alpha(2:end,:)', beta' );                         % Txn
[~,idx] = max(posterior');

% calculate data log-likelihood for forward-backward alg.
data_ll_fb = log(sum(alpha(end,:)));
fprintf('\nFB data log-likelihood = %f\n\n', data_ll_fb);

%%% Mann log-weighted (numerically-stable) forward-backward alg.
eln_alpha = forward_eln( px0, trans_prob, obs_prob, y_obs );
eln_beta = backward_eln( trans_prob, obs_prob, y_obs );
% get the posterior distribution
eln_posterior = elnfb_posterior( eln_alpha(2:end,:)', eln_beta' );
[~,eln_idx] = max(eln_posterior');

% calculate data log-likelihood for ext-log forward-backward alg.
data_ll_elnfb = nansum(nansum(eln_alpha,2));
fprintf('\nExt-Log FB data log-likelihood = %f\n\n', data_ll_elnfb);

%% Liklihood-Weighted Sampling

n = size(px0,1);    % number of states
Ns = 10000;           % number of Monte Carlo sample sequences

% get Ns Monte Carlo sample sequences of length T, and 
% corresponding sequence weights
T = size(y_obs,1);
[ lw_samples, weights ] = lw_sampling( Ns, T, px0, trans_prob, ...
                             obs_prob, y_obs);
                         
% get likelihood-weighted (approximate?) inference posterior
lw_posterior = lw_inference( n, lw_samples, weights );
[~,lw_idx] = max(lw_posterior');

%% Plot Data

% create timeseries
t = linspace(1,size(y_obs,1),5000)';

% create empty continuous state vectors
state     = zeros(size(t,1),1);
eln_state = zeros(size(t,1),1);
lw_state  = zeros(size(t,1),1);

% create plot-friendly continuous states
for i=1:length(t)
    state(i,1)     = idx(floor(t(i)));
    eln_state(i,1) = eln_idx(floor(t(i)));
    lw_state(i,1)  = lw_idx(floor(t(i)));
end

figure(1)
plot(idx,'.','MarkerSize',10); hold on;
plot(eln_idx,'o','MarkerSize',10);
plot(lw_idx,'diamond','MarkerSize',7);
plot(t,state);
plot(t,eln_state,'--');
% plot(t,lw_state,'--');
xlim([0,size(y_obs,1)+1])
ylim([0.5,4.5]); yticks([1 2 3 4 5])
title('HMM - Long Sequence','Interpreter','latex');
xlabel('timestep, $k$','Interpreter','latex');
ylabel('state, $x_k$','Interpreter','latex');
% for problem 2
% hdl = legend('forward-backward','likelihood-weighted inference', ...
%     'fwd-bck state trace', 'likelihood-weighted trace');
% for problem 3
% hdl = legend('forward-backward','ext-log forward-backward', ...
%     'fwd-bck state trace', 'ext-log fwd-bck state trace');
hdl = legend('forward-backward','ext-log forward-backward', ...
    'likelihood-weighted inference','fwd-bck state trace', ...
    'ext-log fwd-bck state trace');
set(hdl,'Interpreter','latex','Location','Northwest')

