function [ data_ll ] = get_data_log_likelihood_eln( eln_alpha )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n = size(eln_alpha,1);
data_ll = 0;

% sum over alpha(X_T)
for i=1:n
    data_ll = data_ll + eexp(eln_alpha(i,end));
end

% return the log of the computed sum;
data_ll = eln(data_ll);

end

