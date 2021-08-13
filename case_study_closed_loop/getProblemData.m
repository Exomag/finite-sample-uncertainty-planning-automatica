function problemData = getProblemData(params)
    % Obstacles' process noise realizations
    problemData.obstacleNoiseSamples = sampleObstacleProcessNoise(params.obstacles, params.Ns, params.N);
    
    % Obstacles' state realizations
    problemData.obstacleStateSamples = sampleObstacleStates(params.obstacles, ...
        problemData.obstacleNoiseSamples, params.N);
    
    % Obstacles' faces realizations
    problemData.obstacleFaceCoeffSamples = sampleObstacleFaceCoeffs(params.obstacles, ...
        problemData.obstacleStateSamples, params.N);

    % Obstacles' faces coefficients moments
    nObstacles = length(params.obstacles);
    obstacleFaceCoeffMean = cell(nObstacles, params.N+1);
    obstacleFaceCoeffVar = cell(nObstacles, params.N+1);
    nConstraints = 0;
    for iObs = 1 : nObstacles
        currObstacle = params.obstacles{iObs};
        nConstraints = nConstraints + currObstacle.numFaces;
        for iHor = 1 : params.N+1
            obstacleFaceCoeffMean{iObs,iHor} = cell(currObstacle.numFaces, 1);
            obstacleFaceCoeffVar{iObs,iHor} = cell(currObstacle.numFaces, 1);
            for iFace = 1 : currObstacle.numFaces
                obstacleFaceCoeffMean{iObs,iHor}{iFace} = ...
                    mean(problemData.obstacleFaceCoeffSamples{1}{iFace}(:,:,iHor), 2);
                obstacleFaceCoeffVar{iObs,iHor}{iFace} = ...
                    cov(problemData.obstacleFaceCoeffSamples{1}{iFace}(:,:,iHor)');
            end
        end
    end
    problemData.nConstraints = nConstraints;
    
    % Constraints' coefficients moments
    constrCoefsSampleMean = cell(params.N+1, nConstraints);
    constrCoefsSampleVar = cell(params.N+1, nConstraints);
    iConstr = 1;
    constrInd2ObsFace = zeros(nConstraints, 2);
    for iObs = 1 : nObstacles
        currObstacle = params.obstacles{iObs};
        for iFace = 1 : currObstacle.numFaces
            for iHor = 1 : params.N+1
                constrCoefsSampleMean{iHor,iConstr} = obstacleFaceCoeffMean{iObs,iHor}{iFace};
                constrCoefsSampleVar{iHor,iConstr} = obstacleFaceCoeffVar{iObs,iHor}{iFace};
            end
            constrInd2ObsFace(iConstr,:) = [iObs iFace];
            iConstr = iConstr + 1;
        end
    end
    problemData.constrInd2ObsFace = constrInd2ObsFace;
    problemData.constrCoefsSampleMean = constrCoefsSampleMean;
    problemData.constrCoefsSampleVar = constrCoefsSampleVar;

    % Concentration bounds
    r1 = zeros(params.N+1, nConstraints);
    r2 = zeros(params.N+1, nConstraints);
    beta_r1 = params.beta / 2 / nObstacles / params.N;
    beta_r2 = params.beta / 2 / nObstacles / params.N;
    for iHor = 1 : params.N+1
        for iConstr = 1 : nConstraints
            only1D = true;
            for i = 1 : size(constrCoefsSampleVar{2}, 1)
                if constrCoefsSampleVar{2} ~= 0
                    only1D = false;
                    break;
                end
            end
            if ~only1D
                lambda = eig(inv(constrCoefsSampleVar{iHor,iConstr}));
                Tval = (3*(params.Ns-1)) / (params.Ns-3) * finv(1-beta_r1, 3, params.Ns-3);
                r1(iHor,iConstr) = real(sqrt(Tval / params.Ns / min(lambda)));
            else
                r1(iHor,iConstr) = tinv(1-beta_r2/2, params.Ns-1) * ...
                    sqrt(constrCoefsSampleVar{iHor,iConstr}(end,end)) / sqrt(params.Ns);
            end
            r2(iHor,iConstr) = max(abs(1 - (params.Ns-1) / chi2inv(beta_r2/2,   params.Ns-1)), ...
                                   abs(1 - (params.Ns-1) / chi2inv(1-beta_r2/2, params.Ns-1)));
        end
    end
    problemData.r1 = r1;
    problemData.r2 = r2;

    % Risk allocation
    problemData.epsilonAlloc = params.epsilon / params.N / nObstacles * ones(params.N, nConstraints);
end
