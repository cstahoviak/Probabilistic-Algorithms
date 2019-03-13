function [ eln_xi, xi ] = elnxi( eln_alpha, eln_beta, trans_prob, obs_prob, y_obs )
%ELNXI - Computes the log of xi(i,j,k)
%   xi(i,j,k) is the probability of being in state x_i and timestep k, and
%   state x_j at timestep k+1 given a model lambda and observation sequence
%   O (in Rabiner notation)

% NOTE: in accordance with Rabiner/Mann notation, the trans_prob matrix
% must be a row-stochastic matrix  

n = size(trans_prob,1);
T = size(y_obs,1);

% initialization
eln_xi = zeros(n,n,T);
xi = zeros(n,n,T);

for k=1:T
    normalizer = NaN;
    for i=1:n
        for j=1:n
            eln_xi(i,j,k) = elnprod( eln_alpha(i,k), ...
                elnprod( eln(trans_prob(i,j)), ...
                elnprod( eln(obs_prob(y_obs(k),j)), eln_beta(j,k+1) )));
            normalizer = elnsum( normalizer, eln_xi(i,j,k) );
        end
    end
    
    for i=1:n
        for j=1:n
            eln_xi(i,j,k) = elnprod( eln_xi(i,j,k), -normalizer);
            xi(i,j,k) = eexp(eln_xi(i,j,k));
        end
    end
end

end

