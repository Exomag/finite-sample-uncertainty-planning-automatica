function problemData = getProblemData(params)

    %% Wall coefficients realizations
    problemData.wallCoefsSamples = sampleWallCoefs(params.wallCoefsMean, params.wallCoefsVar, params.Ns, ...
        params.N);

    %% Wall coefficients moments
    wallCoefsSampleMean = cell(params.N, params.nWalls);
    wallCoefsSampleVar = cell(params.N, params.nWalls);
    for i = 1 : params.N
        for j = 1 : params.nWalls
            wallCoefsSampleMean{i,j} = 1 / params.Ns * sum(problemData.wallCoefsSamples{i,j})';
            wallCoefsSampleVar{i,j} = cov(problemData.wallCoefsSamples{i,j});
        end
    end
    problemData.wallCoefsSampleMean = wallCoefsSampleMean;
    problemData.wallCoefsSampleVar = wallCoefsSampleVar;

    %% Concentration bounds
    r1 = zeros(params.N, params.nWalls);
    r2 = zeros(params.N, params.nWalls);
    beta_r1 = params.beta / 2 / params.N;
    beta_r2 = params.beta / 2 / params.N;
    for i = 1 : params.N
        for j = 1 : params.nWalls
            lambda = eig(inv(wallCoefsSampleVar{i,j}));
            Tval = (3*(params.Ns-1)) / (params.Ns-3) * finv(1-beta_r1, 3, params.Ns-3);
            r1(i,j) = real(sqrt(Tval / params.Ns / min(lambda)));
            
            r2(i,j) = max( abs(1 - (params.Ns-1) / chi2inv(beta_r2/2,   params.Ns-1)), ...
                           abs(1 - (params.Ns-1) / chi2inv(1-beta_r2/2, params.Ns-1)) ...
                         );
        end
    end
    problemData.r1 = r1;
    problemData.r2 = r2;

    %% Risk allocation
    problemData.epsilonAlloc = params.epsilon / params.N * ones(params.N, params.nWalls);
    
    %% Initial conditions
    problemData.xInit = 1;
    problemData.yInit = 1;

    %% Target
    problemData.xTarget = 8;
    problemData.yTarget = 7;
end
