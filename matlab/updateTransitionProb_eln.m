function [ trans_prob_hat ] = updateTransitionProb_eln( eln_gamma, eln_xi)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

n = size(eln_gamma,1);
T = size(eln_gamma,2) - 1;
D = size(eln_gamma,3);

trans_prob_hat = zeros(n,n);

% get transition probability estimate
for i=1:n
    for j=1:n
        numerator = NaN;
        denominator = NaN;

        for d=1:D
            for k=2:T+1
%                 disp([d, k, xi(i,j,k-1,d)])
                numerator = elnsum(numerator, eln_xi(i,j,k-1,d));
                denominator = elnsum(denominator, eln_gamma(i,k-1,d));
            end % end k
        end % end d

        trans_prob_hat(i,j) = eexp(elnprod(numerator, -denominator));
    end % end j  
end % end i

end