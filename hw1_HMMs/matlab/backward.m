function [ beta, beta2 ] = backward( trans_prob, obs_prob, y_obs )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n = size(trans_prob,1);
T = size(y_obs,1);

% initialization
beta = zeros(n,T);
beta(:,end) = ones(n,1);

beta2 = zeros(n,T);
beta2(:,end) = ones(n,1);

% backward pass
for k=(T-1):-1:1
    % matrix math version... not working
%     beta2(:,k) = beta2(:,k+1)'*trans_prob*diag(obs_prob(y_obs(k+1),:));
    beta2(:,k) = obs_prob(y_obs(k+1),:)*trans_prob*diag(beta2(:,k+1));
    
    for i=1:n
        for j=1:n
            beta(i,k) = beta(i,k) + ( beta(j,k+1) * ...
                trans_prob(j,i) * obs_prob(y_obs(k+1),j) );
        end
    end
    % normalize (do NOT normalize beta values!)
%     beta(:,k) = beta(:,k)./sum(beta(:,k));
end

% return beta with dimensions Txn
beta = beta';
beta2 = beta2';

end