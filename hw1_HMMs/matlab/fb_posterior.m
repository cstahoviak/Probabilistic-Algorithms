function [ posterior ] = fb_posterior( alpha, beta )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if size(alpha,2) ~= size(beta,2)
    error('alpha, beta size mismatch')
    
else
    n = size(alpha,1);
    T = size(alpha,2);
    posterior = zeros(n,T);

    for k=1:T
        % sanity check
        fprintf('sum(alpha(:,%d).*beta(:,%d)) = %e\n', k, ...
            k, sum(alpha(:,k).*beta(:,k)))
        
        posterior(:,k) = (alpha(:,k).*beta(:,k)) / ...
            sum( alpha(:,k).*beta(:,k) );
        
        % normalize
        posterior(:,k) = posterior(:,k)./sum(posterior(:,k));
    end
    
    % return a Txn matrix
    posterior = posterior';

end

