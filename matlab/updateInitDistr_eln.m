function [ init_distr ] = updateInitDistr_eln( eln_gamma )
%UNTITLED Summary of this function goes here
%   Copmute pi_i, the estimated probability of starting in state x_i

n = size(eln_gamma,1);
init_distr = zeros(n,1);

for i=1:n
    init_distr(i,1) = eexp(eln_gamma(i,1));
end

end

