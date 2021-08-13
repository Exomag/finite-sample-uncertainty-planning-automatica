function faceCoeffSamples = sampleObstacleFaceCoeffs(obstacles, obstacleStateSamples, N)
    % Sample obstacles' states
    nObstacles = length(obstacles);
    faceCoeffSamples = cell(nObstacles, 1);
    for iObs = 1 : nObstacles
        currObstacle = obstacles{iObs};
        Ns = size(obstacleStateSamples{iObs}, 2);
        faceCoeffSamples{iObs} = cell(currObstacle.numFaces, 1);
        for iFace = 1 : currObstacle.numFaces
            a = currObstacle.aCoeff(iFace,:)';
            faceCoeffSamples{iObs}{iFace} = zeros(length(a)+1, Ns, N+1);
            for iHor = 1 : N+1
                b = currObstacle.cCoeff(iFace,:) * obstacleStateSamples{iObs}(:,:,iHor) + ...
                        currObstacle.dCoeff(iFace);
                faceCoeffSamples{iObs}{iFace}(:,:,iHor) = [repmat(a, 1, Ns); b];
            end
        end
    end
end
