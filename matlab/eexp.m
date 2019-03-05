function [ out ] = eexp( x )
% EEXP - Extended Exponential function
%   Computes the extended exponential as defined by Mann, 2006

% if x == 'LOGZERO'
if isnan(x)
    out = zeros(size(x,1),size(x,2));
else
    out = exp(x);
end

end