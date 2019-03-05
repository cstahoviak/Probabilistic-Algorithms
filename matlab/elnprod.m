function [ prod ] = elnprod( eln_x, eln_y )
% ELNSUM - Extended Logartithm Product function
%   Computes the extended logarithm of the product of x and y given as
%   given as inputs the extended logarithm of x and y, as defined by Mann,
%   2006

% if strcmp(eln(x),'LOGZERO') || strcmp(eln(y),'LOGZERO')
if isnan(eln_x) || isnan(eln_y)
%     prod = 'LOGZERO';
    prod = NaN;
else
    prod = eln_x + eln_y;
end

end