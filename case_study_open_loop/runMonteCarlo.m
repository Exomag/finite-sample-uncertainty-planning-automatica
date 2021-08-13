function pViol = runMonteCarlo(params, problemData, solution)
    %% Parameters
    N_MonteCarlo = 1e5;

    %% Initialization
    Nviol = 0;
    rng(0);

    %% Realizations
    wallCoefsSamples = sampleWallCoefs(params.wallCoefsMean, params.wallCoefsVar, N_MonteCarlo, params.N);

    %% Monte Carlo
    for k = 1 : N_MonteCarlo
        foundCollision = false;
        for i = 1 : params.N
            x_tilde = [solution.xVal(1,i+1) ; solution.xVal(2,i+1) ; 1];
            for j = 1 : params.nWalls
                if wallCoefsSamples{i,j}(k,1) * x_tilde(1) + wallCoefsSamples{i,j}(k,2) * x_tilde(2) + ...
                   wallCoefsSamples{i,j}(k,3) + params.bigM * solution.zVal(j,i) < 0
                    foundCollision = true;
                    Nviol = Nviol + 1;
                    break;
                end
            end
            if foundCollision
                break
            end
        end
    end

    %% Empirical violation
    pViol = Nviol / N_MonteCarlo;
end
