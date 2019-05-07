function [ model, inlier_idx, scores ] = mlesac( data, n, p, t, sigma_vr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

radar_doppler = data(:,1);
radar_azimuth = data(:,2);

Ntargets = size(data,1);
maxIterations = (Ntargets-1)*Ntargets/2;

bestScore = -1e6;   % large negative number
bestInliers = [];
bestModel = [];

scores = [];
k = [];

for i=1:Ntargets-1
    for j=i+1:Ntargets
        
        doppler = [radar_doppler(i); radar_doppler(j)];
        azimuth = [radar_azimuth(i); radar_azimuth(j)];
        
        is_valid = is_data_valid( [doppler, azimuth] );
        if is_valid
            % valid pair of targets 'sampled' - fit the model
            model = doppler2BodyFrameVelocities( doppler', ...
                azimuth', 100);
            
            % generate predicted doppler measurements given the model
            doppler_predicted = simulateRadarDoppler2D( model, ...
                radar_azimuth, zeros(Ntargets,1), zeros(Ntargets,1));
            
            % evaluate the data log-likelihood given the model
            eps_sq = (doppler_predicted - radar_doppler).^2;
            score = -1/(2*sigma_vr^2)*sum(eps_sq);
            
            if score > bestScore
                % this model better explains the data
                bestScore = score;
                bestInliers = (sqrt(eps_sq) < t);
                bestModel = model;
                
                scores = [scores; score];
                
                % evaluate stopping criteria - not yet used
                Ninliers = sum(bestInliers);
                w = Ninliers/Ntargets;
%                 k1 = log(1-0.95)*log(1-w^2);
%                 k2 = log(1-0.75)*log(1-w^2);
%                 k3 = log(1-0.50)*log(1-w^2);
            end
            
        else
            % do nothing - cannot derive a valid model from two targets in
            % the same azimuth bin
        end
        
    end
end

model = bestModel;
inlier_idx = bestInliers;

% % options = optimoptions('lsqnonlin','Display','none');
% f = @(X) ols_func(X, radar_doppler(bestInliers), radar_azimuth(bestInliers));
% 
% % solve OLS problem on inlier set
% model = lsqnonlin(f,bestModel);
% inlier_idx = bestInliers;

end

function [ is_valid ] = is_data_valid( data )

radar_doppler = data(:,1);
radar_azimuth = data(:,2);

numAzimuthBins = getNumAngleBins( radar_azimuth' );

if numAzimuthBins > 1
    is_valid = true;
else
    is_valid = false;
end


end

