function [ alpha, alpha2 ] = forward( px0, trans_prob, obs_prob, y_obs )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n = size(px0,1);
T = size(y_obs,1);

% initialization - alpha(x0) = prior
alpha = zeros(n,T+1);
alpha(:,1) = px0;

alpha2 = zeros(n,T+1);
alpha2(:,1) = px0;

% forward pass
for k=1:T
    % matrix math version... works!
    alpha2(:,k+1) = alpha2(:,k)'*trans_prob'*diag(obs_prob(y_obs(k),:));
    for i=1:n
        for j=1:n
            alpha(i,k+1) = alpha(i,k+1) + ( alpha(j,k) * ...
                trans_prob(i,j) * obs_prob(y_obs(k),i) );
            
%             fprintf('alpha(%d,%d) = \t\t%f\n', j, k, alpha(j,k))
%             fprintf('trans_prob(%d,%d) = \t%f\n', i, j, trans_prob(i,j))
%             fprintf('obs_prob(y_obs(%d),%d) = \t%f\n\n', k, i, obs_prob(y_obs(k),i))
        end
    end
    % normalize (do NOT normalize alpha values!)
%     alpha(:,k+1) = alpha(:,k+1)./sum(alpha(:,k+1));
end

% return alpha with dmensions Txn
alpha = alpha';
alpha2 = alpha2';

end
