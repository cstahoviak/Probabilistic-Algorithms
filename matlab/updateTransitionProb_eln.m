function [ trans_prob_hat ] = updateTransitionProb_eln( eln_gamma, eln_xi )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% NOTE: in accordance with Rabiner/Mann notation, the estimated trans_prob
% matrix will be a row-stochastic matrix

n = size(eln_gamma,1);
T = size(eln_gamma,2);

% initilization
trans_prob_hat = zeros(n,n);

for i=1:n
    for j=1:n
        numerator = NaN;
        denominator = NaN;
        
        for k=1:(T-1)
            numerator = elnsum( numerator, eln_xi(i,j,k) );
            denominator = elnsum( denominator, eln_gamma(i,k) );
        end
        trans_prob_hat(i,j) = eexp(elnprod(numerator, -denominator));
    end
end

end

