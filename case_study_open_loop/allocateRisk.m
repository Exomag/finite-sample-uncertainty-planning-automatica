function finalSolution = allocateRisk(params, problemData, method, solution)
    %% Parameters
    costConvergenceRatio = 1.01;
    activeConstraintThreshold = 1e-3;
    
    %% Initialization
    latestSolution = solution;
    xVal = solution.xVal;
    zVal = solution.zVal;
    totalTime = solution.diagnostic.solvertime;
    prevObj = +inf;
    Nactive = NaN;

    %% Risk allocation
    while (prevObj >= costConvergenceRatio * latestSolution.objVal) && (Nactive ~= 0) && (Nactive ~= params.N)
        % Initialization
        activeConstraints = false(params.N, params.nWalls);
        newRiskSum = 0;
        
        % Check constraints
        for i = 1 : params.N
            x_tilde = [xVal(1,i+1) ; xVal(2,i+1) ; 1];
            for j = 1 : params.nWalls
                if zVal(j,i) == 0
                    if strcmp(method, 'ExactMoments')
                        constraintVal = norminv(1-problemData.epsilonAlloc(i,j)) * ...
                            sqrt(x_tilde' * params.wallCoefsVar{j} * x_tilde) - ...
                            params.wallCoefsMean{j}' * x_tilde;
                    elseif strcmp(method, 'MomentsRobust')
                        constraintVal = norminv(1-problemData.epsilonAlloc(i,j)) * ...
                            sqrt(1+problemData.r2(i,j)) * ...
                            sqrt(x_tilde' * problemData.wallCoefsSampleVar{i,j} * x_tilde) + ...
                            problemData.r1(i,j) * sqrt(x_tilde' * x_tilde) - ...
                            problemData.wallCoefsSampleMean{i,j}' * x_tilde;
                    else
                        error('Solve method not recognized.')
                    end
                    if abs(constraintVal) <= activeConstraintThreshold
                        activeConstraints(i,j) = true;
                        Nactive = Nactive + 1;
                    else
                        if strcmp(method, 'ExactMoments')
                            cdfArg = params.wallCoefsMean{j}' * x_tilde / ...
                                sqrt(x_tilde' * params.wallCoefsVar{j} * x_tilde);
                        elseif strcmp(method, 'MomentsRobust')
                            cdfArg = (problemData.wallCoefsSampleMean{i,j}' * x_tilde - ...
                            problemData.r1(i,j) * sqrt(x_tilde' * x_tilde)) / ...
                            sqrt(x_tilde' * problemData.wallCoefsSampleVar{i,j} * x_tilde) / ...
                            sqrt(1+problemData.r2(i,j));
                        else
                            error('Solve method not recognized.')
                        end
                        newEpsilon = 0.5*problemData.epsilonAlloc(i,j) + 0.5*(1-normcdf(cdfArg));
                        newRiskSum = newRiskSum - (newEpsilon - problemData.epsilonAlloc(i,j));
                        problemData.epsilonAlloc(i,j) = newEpsilon;
                    end
                end
            end
        end
        
        % Reallocate risk
        Nactive = sum(nonzeros(activeConstraints));
        delta = newRiskSum / Nactive;
        for i = 1 : params.N
            for j = 1 : params.nWalls
                if activeConstraints(i,j)
                    problemData.epsilonAlloc(i,j) = problemData.epsilonAlloc(i,j) + delta;
                end
            end
        end
        
        % New solve
        prevObj = latestSolution.objVal;
        latestSolution = solveProblem(params, problemData, method, zVal);
        xVal = latestSolution.xVal;
        totalTime = totalTime + latestSolution.diagnostic.solvertime;
    end
    
    %% Final solution
    finalSolution = latestSolution;
    finalSolution.diagnostic.solvertime = totalTime;
end
