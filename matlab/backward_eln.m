function [ eln_beta ] = backward_eln( trans_prob, obs_prob, y_obs )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% NOTE: in accordance with Rabiner/Mann notation, the trans_prob matrix
% must be a row-stochastic matrix  

n = size(trans_prob,1);
T = size(y_obs,1);

% initialization
eln_beta = zeros(n,T+1);

for k=T:-1:1
    for i=1:n
        logbeta = NaN;
        for j=1:n 
            logbeta = elnsum( logbeta, ...
                elnprod( eln(trans_prob(i,j)), ...
                elnprod( eln(obs_prob(y_obs(k),j)), eln_beta(j,k+1) )));
        end
        eln_beta(i,k) = logbeta;
    end
end

return;


% initialization
eln_beta = zeros(n,T);

for k=(T-1):-1:1
    for i=1:n
        logbeta = NaN;
        for j=1:n 
            logbeta = elnsum( logbeta, ...
                elnprod( eln(trans_prob(i,j)), ...
                elnprod( eln(obs_prob(y_obs(k+1),j)), eln_beta(j,k+1) )));
        end
        eln_beta(i,k) = logbeta;
    end
end

end