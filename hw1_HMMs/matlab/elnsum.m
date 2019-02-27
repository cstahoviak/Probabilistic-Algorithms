function [ sum ] = elnsum( eln_x, eln_y )
% ELNSUM - Extended Logartithm Sum function
%   Computes the extended logarithm of the sum of x and y given as inputs 
%   the extended logarithm of x and y, as defined by Mann, 2006

% if strcmp(eln(x),'LOGZERO') || strcmp(eln(y),'LOGZERO')
if isnan(eln_x) || isnan(eln_y)
%     if strcmp(eln(x),'LOGZERO')
    if isnan(eln_x)
        sum = eln_y;
    else
        sum = eln_x;
    end
else
    if eln_x > eln_y
        sum = eln_x + eln(1 + exp(eln_y-eln_x));
    else
        sum = eln_y + eln(1 + exp(eln_x-eln_y));
    end
end

end