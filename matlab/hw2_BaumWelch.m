%% Header

% Filename:     hw2_BaumWelch.m
% Author:       Carl Stahoviak
% Date Created: 03/02/2019

clc;
clear;
close ALL;

%% Load Data

% load HMM CPT Parameters
load('nominal_hmm_params.mat')

%% Initialize CPT distributions for Baum-Welch

% load observation set for Baum-Welch implementation (y_obs)
% NOTE: each column of y_obs is an independent observation sequence
load('nominal_hmm_multi_logs.mat')

n = size(pxk_xkm1,1);   % number of states
T = size(y_obs,1);      % number of observations per data log
D = size(y_obs,2);      % numnber of unique data logs
p = size(pyk_xk,1);     % number of unique emmision symbols

% number of unique data logs to use in E-step (Berkeley notation)
N_logs = 10;
log_sz = [10, 20, 50, 100];

% number of times M-step will be done
N_mstep = 30;

% total data log-likelihood per EM iteration
data_ll_total = zeros(N_mstep,size(log_sz,2));

% total data log-likelihood per EM iteration
data_ll_mean = zeros(N_mstep,size(log_sz,2));

% scatter plot size
sz = 5;


%% Baum-Welch - Attempt 2

% In this attempt, the e-step (resulting in the calculation of {elnalpha,
% elnbeta, elngamma}) will be run for all D data logs, and a single m-step
% will be done at the end of this. So rather than performing D number of 
% e-step/m-step iterations (which failed at the second iteration), D number
% of e-step iterations will be performed and then a single m-step will be
% done to update the CPT parameter estimates.

figure(1)
ax1 = gca;
title('Baum-Welch EM Parameter Estimation','Interpreter','latex');
xlabel('EM iteration, $s$','Interpreter','latex');
ylabel('data log-likelihood','Interpreter','latex');
hold on;

for i=1:size(log_sz,2)
    
    N_logs = log_sz(i);
    
    % initialize CPTs
    trans_prob = zeros(n,n,N_mstep+1);      % row-stochastic
    obs_prob   = zeros(p,n,N_mstep+1);      % column-stochastic (is this correct?)
    init_distr = zeros(n,N_mstep+1);        % column-stochastic

    % initialize Conditional Probability Tables (CPTs)
    type = 'informed';
    [ init_distr(:,1), trans_prob(:,:,1), obs_prob(:,:,1) ] = initCPTs( ...
        pxk_xkm1', pyk_xk, px0, type);
    
    % data log-likelihood (for each data log per EM iteration)
    data_ll = zeros(N_logs,N_mstep);

    for s=1:N_mstep  % number of times M-step will be done

        fprintf('EM iteration: %d\n', s);

        % init variables
        eln_alpha = zeros(n,T+1,N_logs);    % T+1 columns because alpha(0) is computed
        eln_beta  = zeros(n,T+1,N_logs);    % T+1 columns because beta(0) is computed
        eln_gamma = zeros(n,T+1,N_logs);    % log-posterior
        gamma     = zeros(n,T+1,N_logs);    % true posterior
        eln_xi    = zeros(n,n,T,N_logs);    % log of probability xi(i,j,k)
        xi        = zeros(n,n,T,N_logs);

        for d=1:N_logs  % iterate over D number of data logs
            obs_seq = y_obs(:,d);

            % do E-step for each data log
            [eln_alpha(:,:,d), eln_beta(:,:,d), eln_gamma(:,:,d), ...
                gamma(:,:,d), eln_xi(:,:,:,d), xi(:,:,:,d)] = baumwelch_Estep( ...
                init_distr(:,s), trans_prob(:,:,s), obs_prob(:,:,s), obs_seq);

            % calculate data log-likelihood for the e-step given the
            % current CPT parameter estimates
            data_ll(d,s) = get_data_log_likelihood_eln( eln_alpha(:,:,d) );
        end

        % do M-step given D iterations of the E-step
        [init_distr(:,s+1), trans_prob(:,:,s+1), obs_prob(:,:,s+1) ] = ...
            baumwelch_Mstep( p, eln_gamma, gamma, eln_xi, xi, y_obs);

        % plot the data log-likelihood of each data log at EM iteration s
        scatter(ax1,s*ones(N_logs,1),data_ll(:,s),sz,'filled');

        % compute the total data log-likelihood per EM iteration
        data_ll_total(s,i) = sum(data_ll(:,s));
        data_ll_mean(s,i) = mean(data_ll(:,s));

    end

end

fprintf('\nFinal (sum) data log-likelihood =\n')
disp([log_sz; data_ll_total(end,:)]);
fprintf('Final (avg) data log-likelihood =\n\n')
disp([log_sz; data_ll_mean(end,:)]);

fprintf('Estimated Initial Distribution\n')
disp(init_distr(:,end));
disp(sum(init_distr(:,end),1));

fprintf('True State Transition CPT\n')
disp(pxk_xkm1);
fprintf('%d) Estimated State Transition CPT\n', s)
disp(trans_prob(:,:,end)');
disp(1-sum(trans_prob(:,:,end)',1));

fprintf('True Emission CPT\n')
disp(pyk_xk);
fprintf('Estimated Emission CPT\n')
disp(obs_prob(:,:,end));
disp(1-sum(obs_prob(:,:,end),1));

% fprintf('State Transition CPT Residual:\n')
% disp(pxk_xkm1 - trans_prob(:,:,end)');
% 
% fprintf('Emission CPT Residual:\n')
% disp(pyk_xk - obs_prob(:,:,end));

figure(2)
for i=1:size(log_sz,2)
    plot(data_ll_mean(:,i));
    title('Baum-Welch EM Parameter Estimation','Interpreter','latex');
    xlabel('EM iteration, $s$','Interpreter','latex');
    ylabel('average data log-likelihood','Interpreter','latex');
    hold on;
    hdl = legend();
end
set(hdl,'Interpreter','latex','Location','best');

return;

%% Baum-Welch - Attempt 1 - INCORRECT

% In this implementation, I attempted to run the e-step immediately
% followed by the m-step for a single observation sequence. This way, an
% updater CPT parameter set {trans_prob, obs_prob, and init_distr} was
% estimated at the end of a single observation sequence (or data log). This
% method is incorrect... for reasons explained by Nisar and the Berkeley
% paper.

% initialize
eln_xi = zeros(n,n,T-1,D);      % log of probability xi(i,j,k)
xi = zeros(n,n,T-1,D);

for i=1:D
    obs_seq = y_obs(:,i);

    %%% Mann log-weighted (numerically-stable) forward-backward alg.
    eln_alpha = forward_eln( init_distr(:,i), ...
        trans_prob(:,:,i), obs_prob(:,:,i), obs_seq );
    eln_beta = backward_eln( trans_prob(:,:,i), ...
        obs_prob(:,:,i), obs_seq );
    
    fprintf('eln_alpha(%d) =\n',i);
    disp([(0:T)' eln_alpha'])
    fprintf('eln_beta(%d) =\n', i);
    disp([(1:T)' eln_beta'])

    % get the log-posterior and posterior distributions (exact inference)
    [ eln_gamma, gamma ] = posterior_elnfb( eln_alpha(:,2:end), eln_beta );
    
    fprintf('\neln_gamma(%d) =\n', i);
    disp([(1:T)' eln_gamma'])
    
    % get log of probability xi(i,j,k)
    [ eln_xi(:,:,:,i), xi(:,:,:,i) ] = elnxi( eln_alpha(:,2:end), ...
        eln_beta, trans_prob(:,:,i), obs_prob(:,:,i), obs_seq );

    % update state transition probability, pxk_xkm1
    % NOTE: trans_prob is a row-stochastic matrix (Rabiner notation)
    [ trans_prob(:,:,i+1) ] = updateTransitionProb_eln( eln_gamma, ...
        eln_xi(:,:,:,i) );
    
    fprintf('%d) True State Transition CPT\n', i)
    disp(pxk_xkm1);
    fprintf('%d) Estimated State Transition CPT\n', i)
    disp(trans_prob(:,:,i+1)');
    disp(sum(trans_prob(:,:,i+1)',1));
    
    % update emission probability, pyk_xk
    [ obs_prob(:,:,i+1) ] = updateEmissionProb_eln( p, ... 
        eln_gamma, eln_xi(:,:,:,i), obs_seq );
    
    fprintf('%d) True Emission CPT\n', i)
    disp(pyk_xk);
    fprintf('%d) Estimated Emission CPT\n', i)
    disp(obs_prob(:,:,i+1));
    disp(sum(obs_prob(:,:,i+1),1));
    
    % update initial distribution, px0
    init_distr(:,i+1) = updateInitDistr_eln( eln_gamma );
    
    fprintf('%d) Estimated Initial Distribution\n', i)
    disp(init_distr(:,i+1));
    disp(sum(init_distr(:,i+1),1));
    
end

% fprintf('State Transition CPT Residual:\n')
% disp(pxk_xkm1 - trans_prob(:,:,end));
% 
% fprintf('Emission CPT Residual:\n')
% disp(pyk_xk - obs_prob(:,:,end));

