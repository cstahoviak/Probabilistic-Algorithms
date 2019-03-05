function [ init_distr ] = estimateInitDistr( eln_gamma )
%UNTITLED Summary of this function goes here
%   Copmute pi_i, the estimated probability of starting in state x_i

init_distr = eexp(eln_gamma(:,1));

end

