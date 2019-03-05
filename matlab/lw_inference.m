function [ posterior ] = lw_inference( n, lw_samples, weights )
%LW_INFERENCE Likelihood-weighted approximate inference
%   Detailed explanation goes here

T = size(lw_samples,2);
posterior = zeros(n,T);

for k=1:T
    for i=1:n
        % get index location of where x_k = x_i across all sequences
        idx = (lw_samples(:,k) == i);
        
        % the posterior for each state i at timestep k is the weighted sum
        % of the realizations (indicator function) of state x_i at
        % timestep k across all MC sample sequences
        posterior(i,k) = sum(weights(idx))/sum(weights);
    end
end

posterior = posterior';

end

