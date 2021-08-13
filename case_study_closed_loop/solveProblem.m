function solution = solveProblem(params, problemData)

    %% Variables
    xVar = sdpvar(params.nx, params.N+1);
    uVar = sdpvar(params.nu, params.N);
    zVar = binvar(problemData.nConstraints, params.N);

    %% Objective
    objective = 0;
    for i = 2 : params.N+1
        objective = objective + ...
            uVar(:,i-1)' * params.inputWeight * uVar(:,i-1) + ...
            (xVar(:,i) - params.targetState)' * params.targetStateWeight * (xVar(:,i) - params.targetState);
    end

    %% Equality constraints
    constraints = (xVar(:,1) == params.egoInit);
    for i = 1 : params.N
        constraints = [constraints, xVar(:,i+1) == params.A_discrete*xVar(:,i) + params.B_discrete*uVar(:,i)];
    end

    %% Inequality constraints
    for i = 1 : params.N
        % Position
        constraints = [constraints, params.yMin + params.egoWidth/2 <= xVar(3,i+1) ...
            <= params.yMax - params.egoWidth / 2];

        % Acceleration
        constraints = [constraints, params.xAccMin <= uVar(1,i) <= params.xAccMax];
        constraints = [constraints, params.yAccMin <= uVar(2,i) <= params.yAccMax];
        
        % Obstacle constraints
        for j = 1 : problemData.nConstraints
            for k = 1 : 5
                switch k
                    case 1
                        x_tilde = [xVar(1,i+1) ; xVar(3,i+1) ; 1];
                    case 2
                        x_tilde = [xVar(1,i+1)+params.egoLength/2 ; xVar(3,i+1)+params.egoWidth/2 ; 1];
                    case 3
                        x_tilde = [xVar(1,i+1)+params.egoLength/2 ; xVar(3,i+1)-params.egoWidth/2 ; 1];
                    case 4
                        x_tilde = [xVar(1,i+1)-params.egoLength/2 ; xVar(3,i+1)+params.egoWidth/2 ; 1];
                    case 5
                        x_tilde = [xVar(1,i+1)-params.egoLength/2 ; xVar(3,i+1)-params.egoWidth/2 ; 1];
                end
                constraints = [constraints, norminv(1-problemData.epsilonAlloc(i,j)) * ...
                    sqrt(1+problemData.r2(i+1,j)) * ...
                    sqrt(x_tilde' * problemData.constrCoefsSampleVar{i+1,j} * x_tilde) + ...
                    problemData.r1(i+1,j) * norm(x_tilde, 2) <= ...
                    problemData.constrCoefsSampleMean{i+1,j}' * x_tilde + params.bigM * zVar(j,i)];
            end
        end
        for iObs = 1 : length(params.obstacles)
            currObstacle = params.obstacles{iObs};
                constraints = [constraints, sum(zVar(problemData.constrInd2ObsFace(:,1)==iObs,i)) <= ...
                    currObstacle.numFaces-1];
        end
    end

    %% Solve
    options = sdpsettings('verbose', 0, 'solver', 'cplex');
    diagnostic = optimize(constraints, objective, options);
    solution.diagnostic = diagnostic; 

    %% Solution data
    if diagnostic.problem == 0
        solution.tVal = params.Ts * (0 : params.N);
        solution.xVal = value(xVar);
        solution.uVal = value(uVar);
        solution.zVal = value(zVar);
        solution.objVal = value(objective);
    else
        disp(yalmiperror(diagnostic.problem));
    end
    
end