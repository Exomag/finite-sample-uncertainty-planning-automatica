function params = loadParams()
    %% Horizon
    % Length
    params.N = 10;
    
    % Sampling time
    params.Ts = 1;

    %% Dimensions
    % Corridor
    params.corridorWidth = 2;
    
    %% Constraints
    % Position
    params.xMin = 0;
    params.yMin = 0;
    params.xMax = 9;
    params.yMax = 9;
    
    % Velocity
    params.xVelMin = -1;
    params.yVelMin = -1;
    params.xVelMax = 1;
    params.yVelMax = 1;
    
    % Number of walls
    params.nWalls = 2;
    
    % Wall coefficients mean/variance
    params.wallCoefsMean{1} = [-1; 0; 2];
    params.wallCoefsMean{2} = [0; 1; -6];
    params.wallCoefsVar{1} = 0.001 * eye(3);
    params.wallCoefsVar{2} = 0.001 * eye(3);

    % Wall polygon vertices
    params.wallPolygonVertMean = getWallPolygonVert(params.wallCoefsMean, params.yMin, params.xMax);

    %% Dynamics
    % Dimensions
    params.nx = 2;
    params.nu = 2;
    
    % Continuous dynamics
    params.A_continuous = zeros(params.nx, params.nx);
    params.B_continuous = eye(params.nu);
    params.systemContinuous = ss(params.A_continuous, params.B_continuous, ...
        zeros(params.nu,params.nx), zeros(params.nu));
    
    % Discrete dynamics
    params.systemDiscrete = c2d(params.systemContinuous, params.Ts);
    params.A_discrete = params.systemDiscrete.A;
    params.B_discrete = params.systemDiscrete.B;

    %% Cost weights
    params.weightCost = 1;

    %% Chance constraints
    % Safety margin
    params.epsilon = 0.05;
    
    % Certainty
    params.beta = 1e-3;
    
    % Number of samples
    params.Ns = ceil( exp(1)/(exp(1)-1) / params.epsilon * ...
        (log(2^(params.nWalls*params.N)/params.beta) + params.nu*params.N - 1));

    %% Other
    params.bigM = 1000;
end