function [ lw_samples, weights ] = lw_sampling( Ns, px0, trans_prob, obs_prob, y_obs )
%LW_SAMPLING Liklihood Weighted Sampling
%   Detailed explanation goes here

T = size(y_obs,1);

lw_samples = zeros(Ns,T);   % V_E - likelihood weighted sample
weights = ones(Ns,1);       % W_E - likelihood weights

% get intial state from prior distribution, px0
[~,idx] = max(px0);

for s=1:Ns
    for k=1:T
        if k==1
            % draw random sample according to initial (known) state, x0
            lw_samples(s,k) = randsample(4,1,true, ...
                trans_prob(:,idx)); 
        else
            % draw random sample according to state transtion probabilities
            lw_samples(s,k) = randsample(4,1,true, ...
                trans_prob(:,lw_samples(s,k-1))); 
        end
        
        % update sequence weight:
        % weight = weight * P( y_k | Pa(y_k) = x_k )
        % NOTE: parent of each observation y_k is the state x_k
        weights(s,1) = weights(s,1) * ...
            obs_prob(y_obs(k),lw_samples(s,k));
    end
end

end

