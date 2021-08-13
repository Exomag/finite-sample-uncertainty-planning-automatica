function solution = solveProblem(params, problemData, method, zVal)
    %% Variables
    xVar = sdpvar(params.nx, params.N+1);
    uVar = sdpvar(params.nu, params.N);
    if nargin < 4
        noBinaries = false;
        zVar = binvar(params.nWalls, params.N);
    else
        noBinaries = true;
    end

    %% Objective
    target = [problemData.xTarget; problemData.yTarget];
    objective = params.weightCost * (xVar(:,end)-target)' * (xVar(:,end)-target);

    %% Equality constraints
    x_0 = [problemData.xInit; problemData.yInit];
    constraints = (xVar(:,1) == x_0);
    for i = 1 : params.N
        constraints = [constraints, xVar(:,i+1) == params.A_discrete*xVar(:,i) + params.B_discrete*uVar(:,i)];
    end

    %% Inequality constraints
    for i = 1 : params.N
        % Position
        constraints = [constraints, xVar(1,i+1) >= params.xMin];
        constraints = [constraints, xVar(1,i+1) <= params.xMax];
        constraints = [constraints, xVar(2,i+1) >= params.yMin];
        constraints = [constraints, xVar(2,i+1) <= params.yMax];

        % Velocity
        constraints = [constraints, uVar(1,i) >= params.xVelMin];
        constraints = [constraints, uVar(1,i) <= params.xVelMax];
        constraints = [constraints, uVar(2,i) >= params.yVelMin];
        constraints = [constraints, uVar(2,i) <= params.yVelMax];
        
        % Walls
        x_tilde = [xVar(1,i+1) ; xVar(2,i+1) ; 1];
        for j = 1 : params.nWalls
            if ~noBinaries
                bigMTerm = params.bigM * zVar(j,i);
            else
                if zVal(j,i) == 1
                    continue;
                end
                bigMTerm = 0;
            end
            if strcmp(method, 'ExactMoments')
                constraints = [constraints, norminv(1-problemData.epsilonAlloc(i,j)) * ...
                    sqrt(x_tilde' * params.wallCoefsVar{j} * x_tilde) <= ...
                    params.wallCoefsMean{j}' * x_tilde + bigMTerm];
            elseif strcmp(method, 'MomentsRobust')
                constraints = [constraints, norminv(1-problemData.epsilonAlloc(i,j)) * ...
                    sqrt(1+problemData.r2(i,j)) * ...
                    sqrt(x_tilde' * problemData.wallCoefsSampleVar{i,j} * x_tilde) + ...
                    problemData.r1(i,j) * sqrt(x_tilde' * x_tilde) <= ...
                    problemData.wallCoefsSampleMean{i,j}' * x_tilde + bigMTerm];
            elseif strcmp(method, 'Scenario')
                constraints = [constraints, problemData.wallCoefsSamples{i,j}(:,1) * x_tilde(1) + ...
                    problemData.wallCoefsSamples{i,j}(:,2) * x_tilde(2) + ...
                    problemData.wallCoefsSamples{i,j}(:,3) + bigMTerm >= 0];
            else
                error('Solve method not recognized.')
            end
        end
        if ~noBinaries
            constraints = [constraints, sum(zVar(:,i)) <= params.nWalls-1];
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
        if ~noBinaries
            solution.zVal = value(zVar);
        else
            solution.zVal = zVal;
        end
        solution.objVal = value(objective);
    else
        disp(yalmiperror(diagnostic.problem));
    end
end
