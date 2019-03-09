function [ obs_prob_hat ] = updateEmissionProb_eln( p, eln_gamma, y_obs )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

n = size(eln_gamma,1);
T = size(eln_gamma,2) - 1;
D = size(y_obs,2);

obs_prob_hat = zeros(p,n);

% get emission probabilty estimate
for j=1:p
    for i=1:n
        numerator = NaN;
        denominator = NaN;

        for d=1:D
            for k=2:T+1
%                 disp(k)
                % symbol j observed at time step k in data log d
                if y_obs(k-1,d) == j
                    numerator = elnsum(numerator, eln_gamma(i,k,d));
                    if j == 15
%                         disp('HERE')
%                         disp([j, k, d])
                    end
                end
                denominator = elnsum(denominator, eln_gamma(i,k,d));
            end % end k
        end % end d

        obs_prob_hat(j,i) = eexp(elnprod(numerator, -denominator));
    end % end i
end % end j

end