%% Header

% Filename:     hw2_mle.m
% Author:       Carl Stahoviak
% Date Created: 03/02/2019

clc;
clear;
close ALL;

%% Load Data

% load HMM CPT Parameters
load('nominal_hmm_params.mat')
disp(pxk_xkm1)

%% Mann Extended-Logarithm Forward-Backward Algorithm

% load observation set for Baum-Welch implementation (y_obs)
% NOTE: each column of y_obs is an independent observation sequence
load('nominal_hmm_multi_logs.mat')

T = size(y_obs,1);      % number of observations
n = size(pxk_xkm1,1);   % number of states
p = size(pyk_xk,1);     % number of unique emmision symbols

% number of unique data logs to use in estimation of CPTs
% N = size(y_obs,2);  % N \in [1,size(y_obs,2)]
N = 2;

% initialize CPTs
eln_xi = zeros(n,n,T-1,N);
trans_prob = zeros(n,n,N+1);    % row-stochastic
obs_prob = zeros(p,n,N+1);      % column-stochastic (is this correct?)

% initialize CPTs (pxk_xkm1 and pyk_xk)
trans_prob(:,:,1) = pxk_xkm1';  % row-stochastic
obs_prob(:,:,1)   = pyk_xk;
init_distr        = (1/n)*ones(n,1);

% try a different initialization...
trans_prob(:,:,1) = [1/3 0   1/3 1/3;
                     1/4 1/4 1/4 1/4;
                     0   1/3 1/3 1/3;
                     1/2 0   0   1/2];

for i=1:N
    obs_seq = y_obs(:,i);

    %%% Mann log-weighted (numerically-stable) forward-backward alg.
    eln_alpha = forward_eln( init_distr, trans_prob(:,:,i), ...
        obs_prob(:,:,i), obs_seq );
    eln_beta = backward_eln( trans_prob(:,:,i), ...
        obs_prob(:,:,i), obs_seq );
    
    fprintf('eln_alpha(%d) =\n',i);
    disp([(0:T)' eln_alpha])
    fprintf('eln_beta(%d) =\n', i);
    disp([(1:T)' eln_beta])

    % get the log-posterior and posterior distributions (exact inference)
    [ eln_gamma, gamma ] = posterior_elnfb( eln_alpha(2:end,:)', eln_beta' );
    
    fprintf('\neln_gamma(%d) =\n', i);
    disp([(1:T)' eln_gamma])
    
    % get log of probability xi(i,j,k)
    [ eln_xi(:,:,:,i) ] = elnxi( eln_alpha(2:end,:)', eln_beta', ...
        trans_prob(:,:,i), obs_prob(:,:,i), obs_seq );

    % update state transition probability, pxk_xkm1
    % NOTE: trans_prob is a row-stochastic matrix (Rabiner notation)
    [ trans_prob(:,:,i+1) ] = estimateTransitionProb( eln_gamma', ...
        eln_xi(:,:,:,i) );
    
    fprintf('\n%d) True State Transition CPT\n', i)
    disp(pxk_xkm1);
    fprintf('%d) Estimated State Transition CPT\n', i)
    disp(trans_prob(:,:,i+1)');
    disp(sum(trans_prob(:,:,i+1)',1));
    
    % update emission probability, pyk_xk
    [ obs_prob(:,:,i+1) ] = estimateEmissionProb( p, ... 
        eln_gamma', eln_xi(:,:,:,i), obs_seq );
    
    fprintf('%d) True Emission CPT\n', i)
    disp(pyk_xk);
    fprintf('%d) Estimated Emission CPT\n', i)
    disp(obs_prob(:,:,i+1));
    disp(sum(obs_prob(:,:,i+1),1));
    
end

fprintf('State Transition CPT Residual:\n')
disp(pxk_xkm1 - trans_prob(:,:,end));

fprintf('Emission CPT Residual:\n')
disp(pyk_xk - obs_prob(:,:,end));

return;


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

