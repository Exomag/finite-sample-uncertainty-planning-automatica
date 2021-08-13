function stateSamples = sampleObstacleStates(obstacles, obstacleNoiseSamples, N)
    % Sample obstacles' states
    nObstacles = length(obstacles);
    stateSamples = cell(nObstacles, 1);
    for iObs = 1 : nObstacles
        currObstacle = obstacles{iObs};
        Ns = size(obstacleNoiseSamples{iObs}, 2);
        stateSamples{iObs} = zeros(currObstacle.nx, Ns, N+1);
        for iHor = 1 : N+1
            if iHor == 1
                stateSamples{iObs}(:,:,iHor) = mvnrnd(currObstacle.initMean, currObstacle.initVar, Ns)';
            else
                w = obstacleNoiseSamples{iObs}(:,:,iHor-1);
                stateSamples{iObs}(:,:,iHor) = currObstacle.E * stateSamples{iObs}(:,:,iHor-1) + currObstacle.F + w;
            end
        end
    end
end
