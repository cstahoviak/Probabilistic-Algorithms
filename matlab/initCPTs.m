function [ init_distr, trans_prob, obs_prob ] = ...
    initCPTs( pxk_xkm1, pyk_xk, px0, type )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% NOTE: in accordance with Rabiner/Mann notation, the trans_prob matrix
% must be a row-stochastic matrix  

n = size(pxk_xkm1,1);   % number of states
p = size(pyk_xk,1);     % number of unique emmision symbols

if strcmp(type, 'truth')
    init_distr = px0;
    trans_prob = pxk_xkm1;          % row-stochastic
    obs_prob   = pyk_xk;
    
elseif strcmp(type, 'informed')
    % "informed" uniform initialization
    % initialize CPTs WITH knowledge of where zeros exist in the tables
    init_distr = (1/n)*ones(n,1);

    trans_prob = [1/3 0   1/3 1/3;
                         1/4 1/4 1/4 1/4;
                         0   1/3 1/3 1/3;
                         1/2 0   0   1/2];

    obs_prob   = [(1/10)*ones(10,1), zeros(10,2), (1/10)*ones(10,1);
                  zeros(4,1), (1/4)*ones(4,1), (1/5)*ones(4,1), zeros(4,1);
                  0,          0,               1/5,             0];
                  
elseif strcmp(type, 'uninformed')
    % "un-informed" uniform initialization
    % initialize CPTs WITHOUT knowledge of where zeros exist in the tables
    init_distr = (1/n)*ones(n,1);
    trans_prob = (1/n)*ones(n);
    obs_prob   = (1/p)*ones(p,n);
    
elseif strcmp(type, 'rand')
    init_distr = rand(n,1);
    init_distr = init_distr./sum(init_distr);   % normalize
    
    trans_prob = rand(n);
    for i=1:n
        % normalize rows
        trans_prob(i,1) = trans_prob(i,1)./sum(trans_prob(i,1));
    end
    
    obs_prob = rand(p,n);
    for i=1:n
        % normalize columns
        obs_prob(:,i) = obs_prob(:,i)./sum(obs_prob(:,i));
    end
    
else
    error('Improperly specified initialization type')
end
   
    
end

