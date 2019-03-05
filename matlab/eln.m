function [ out ] = eln( x )
% ELN - Extended Natural Logarathm function
%   Computes the extended natural logarithm as defined by Mann, 2006

if x == 0
%     out = 'LOGZERO';
    out = NaN;
elseif x > 0
    out = log(x);
else
    error('eln() negative input error')   
end

end