function wallCoefsSamples = sampleWallCoefs(wallCoefsMean, wallCoefsVar, Ns, N)
    
    %% Wall coefficients realizations
    nWalls = length(wallCoefsMean);
    wallCoefsSamples = cell(N, nWalls);
    for i = 1 : N
        for j = 1 : nWalls
            wallCoefsSamples{i,j} = mvnrnd(wallCoefsMean{j}, wallCoefsVar{j}, Ns);
        end
    end
    
end

