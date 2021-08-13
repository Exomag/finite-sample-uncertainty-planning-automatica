function noiseSamples = sampleObstacleProcessNoise(obstacles, Ns, N)
    % Sample obstacles' process noise
    nObstacles = length(obstacles);
    noiseSamples = cell(nObstacles, 1);
    for i = 1 : nObstacles
        currObstacle = obstacles{i};
        noiseSamples{i} = zeros(currObstacle.nx, Ns, N);
        for j = 1 : N
            noiseSamples{i}(:,:,j) = mvnrnd(zeros(1, currObstacle.nx), currObstacle.G, Ns)';
        end
    end
end
