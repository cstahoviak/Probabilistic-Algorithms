function [ eps ] = ols_func( X, radar_doppler, radar_azimuth )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

model = X;
Ntargets = size(radar_doppler,1);

% get predicted 
doppler_predicted = simulateRadarDoppler2D( model, radar_azimuth, ...
    zeros(Ntargets,1), zeros(Ntargets,1));

% define error vector, epsilon
eps = doppler_predicted - radar_doppler;

end

