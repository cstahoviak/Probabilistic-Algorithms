function [ eln_gamma, gamma ] = posterior_elnfb( eln_alpha, eln_beta )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if size(eln_alpha,2) ~= size(eln_beta,2)
    error('eln_alpha, eln_beta size mismatch')
    
else
    n = size(eln_alpha,1);
    T = size(eln_alpha,2);
    
    % initialization
    eln_gamma = zeros(n,T);     % log-posterior
    gamma = zeros(n,T);         % true posterior
    
    for k=1:T
        normalizer = NaN;
        for i=1:n
            eln_gamma(i,k) = elnprod(eln_alpha(i,k),eln_beta(i,k));
            normalizer = elnsum(normalizer,eln_gamma(i,k));
        end
        
        for i=1:n
            eln_gamma(i,k) = elnprod(eln_gamma(i,k),-normalizer);
            gamma(i,k) = eexp(eln_gamma(i,k));
        end
        % sanity check
        fprintf('1 - sum(gamma(:,%d)) = %e\n', k, ...
            1-nansum(eexp(eln_gamma(:,k))))
    end
end

% return a Txn matrix
eln_gamma = eln_gamma';
gamma = gamma';

end

