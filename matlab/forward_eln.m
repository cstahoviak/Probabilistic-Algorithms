function [ eln_alpha ] = forward_eln( px0, trans_prob, obs_prob, y_obs )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% NOTE: in accordance with Rabiner/Mann notation, the trans_prob matrix
% must be a row-stochastic matrix  

n = size(px0,1);
T = size(y_obs,1);
eln_alpha = zeros(n,T+1);

% initialization
for i=1:n    
%     eln_alpha(i,1) = elnprod( eln(px0(i,1)), ...
%         eln(obs_prob(y_obs(1),i)) );
    eln_alpha(i,1) = eln(px0(i,1));
end

for k=2:T+1
    for j=1:n
        logalpha = NaN;
        for i=1:n 
            logalpha = elnsum( logalpha, ...
                elnprod( eln_alpha(i,k-1), eln(trans_prob(i,j)) ));
        end
        eln_alpha(j,k) = elnprod( logalpha, ...
            eln(obs_prob(y_obs(k-1),j)) );
    end
end

end