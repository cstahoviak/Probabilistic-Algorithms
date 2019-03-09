function [ obs_prob_hat] = updateEmissionProb( p, gamma, y_obs )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

n = size(gamma,1);
T = size(gamma,2) - 1;
D = size(gamma,3);

obs_prob_hat = zeros(p,n);

% get emission probabilty estimate
for j=1:p
    for i=1:n
        numerator = 0;
        denominator = 0;

        for d=1:D
            for k=2:T+1
                % symbol j observed at time step k in data log d
                if y_obs(k-1,d) == j
                    numerator = numerator + gamma(i,k,d);
                end
                denominator = denominator + gamma(i,k,d);
            end % end k
        end % end d

        obs_prob_hat(j,i) = numerator/denominator;
    end % end i
end % end j

% normalize columns of emission probability table
% for i=1:n
%     obs_prob_hat(:,i) = obs_prob_hat(:,i)./sum(obs_prob_hat(:,i));
% end

end

