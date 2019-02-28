function [ state, trace ] = getStateTrace( tspan, posterior )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% assumes posterior to be an nxT matrix
[~,idx] = max(posterior);
state = idx;

% create empty continuous state vectors
trace = zeros(size(tspan,1),1);

% create plot-friendly continuous states
for i=1:length(tspan)
    trace(i,1) = idx(floor(tspan(i)));
end

end

