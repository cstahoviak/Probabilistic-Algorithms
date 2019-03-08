function [ trans_prob_hat ] = updateTransitionProb( gamma, xi)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

n = size(gamma,1);
T = size(gamma,2);
D = size(gamma,3);

trans_prob_hat = zeros(n,n);

% get transition probability estimate
for i=1:n
    for j=1:n
        numerator = 0;
        denominator = 0;

        for d=1:D
            for k=2:T
%                 disp([d, k, xi(i,j,k-1,d)])
                numerator = numerator + xi(i,j,k-1,d);
                denominator = denominator + gamma(i,k-1,d);
            end % end k
        end % end d

        trans_prob_hat(i,j) = numerator/denominator;
    end % end j  
end % end i

% normalize rows of transition probability table
% for i=1:n
%     trans_prob_hat(i,:) = trans_prob_hat(i,:)./sum(trans_prob_hat(i,:));
% end
    
end

