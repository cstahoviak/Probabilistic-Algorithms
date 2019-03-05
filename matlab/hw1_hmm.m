%% Header

% Author: Carl Stahoviak
% Date Created:

clc;
clear;
close ALL;

%% Load Data

load('nominal_hmm_params.mat')
trans_prob = pxk_xkm1;  % column-stochastic
obs_prob = pyk_xk;      % column-stochastic

% load('nominal_hmm_short_log.mat')
load('nominal_hmm_long_log.mat')

%% The Forward-Backward Algorithm

%%% standard forward-backward algorithm
[alpha, alpha2] = forward( px0, trans_prob, obs_prob, y_obs );  % Txn
[beta, beta2] = backward( trans_prob, obs_prob, y_obs );        % Txn

% get the posterior distribution (exact inference)
posterior = posterior_fb( alpha(2:end,:)', beta' );             % Txn

% calculate data log-likelihood for forward-backward alg.
data_ll_fb = log(sum(alpha(end,:)));
fprintf('\nFB data log-likelihood = %f\n\n', data_ll_fb);

%% Mann Extended-Logarithm Forward-Backward Algorithm

% use Mann notation
trans_prob = trans_prob';   % row-stochastic (Rabiner/Mann convention)

%%% Mann log-weighted (numerically-stable) forward-backward alg.
eln_alpha = forward_eln( px0, trans_prob, obs_prob, y_obs );
eln_beta = backward_eln( trans_prob, obs_prob, y_obs );

% get the log-posterior distribution (exact inference)
[ eln_gamma, gamma ] = posterior_elnfb( eln_alpha(2:end,:)', eln_beta' );
eln_posterior = gamma;  % true posterior, not log-posterior

% calculate data log-likelihood for ext-log forward-backward alg.
nonNaN_idx = ~isnan(eln_alpha(end,:));
data_ll_elnfb = eln(sum(eexp(eln_alpha(end,nonNaN_idx))));
fprintf('\nExt-Log FB data log-likelihood = %f\n\n', data_ll_elnfb);

%% Liklihood-Weighted Sampling

n = size(px0,1);    % number of states
Ns = 10000;         % number of Monte Carlo sample sequences

% get Ns Monte Carlo sample sequences of length T, and 
% corresponding sequence weights
T = size(y_obs,1);
[ lw_samples, weights ] = lw_sampling( Ns, px0, trans_prob, ...
                            obs_prob, y_obs );
                         
% get likelihood-weighted approximate inference posterior
lw_posterior = lw_inference( n, lw_samples, weights );

%% Plot Data

% create timeseries
tspan = linspace(1,size(y_obs,1),5000)';

% get continuous state trace
% NOTE: assumes posterior to be an nxT matrix
[ state, trace ] = getStateTrace( tspan, posterior' );
[ eln_state, eln_trace ] = getStateTrace( tspan, eln_posterior' );
[ lw_state, lw_trace ] = getStateTrace( tspan, lw_posterior' );

figure(1)
plot(state,'.','MarkerSize',10); hold on;
plot(eln_state,'o','MarkerSize',10);
plot(lw_state,'diamond','MarkerSize',7);
plot(tspan,trace);
plot(tspan,eln_trace,'--');
% plot(tspan,lw_trace,'--');
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

