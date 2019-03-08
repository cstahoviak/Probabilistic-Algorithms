function [ out ] = eexp( x )
% EEXP - Extended Exponential function
%   Computes the extended exponential as defined by Mann, 2006

% if x == 'LOGZERO'
if isnan(x)
    out = 0;
else
    out = exp(x);
end

end