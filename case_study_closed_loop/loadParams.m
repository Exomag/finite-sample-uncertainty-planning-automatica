function params = loadParams()

    %% Horizon
    % Length
    params.N = 25;
    
    % Sampling time
    params.Ts = 0.2;
    
    % Time vector
    params.time = 0 : params.Ts : params.N * params.Ts;

    %% Lanes
    % Lane width
    params.laneWidth = 4;
    
    % Lane center
    params.laneCenter1 = - params.laneWidth / 2;
    params.laneCenter2 = + params.laneWidth / 2;
    
    %% Constraints
    % Position
    params.yMin = 0;
    params.yMax = params.laneWidth;
    
    % Acceleration
    params.xAccMin = -10;
    params.xAccMax = 10;
    params.yAccMin = -2;
    params.yAccMax = 2;

    %% Ego dynamics
    % State-space dimensions
    params.nx = 4;
    params.nu = 2;
    
    % Discrete dynamics
    A_long = [1 params.Ts; 0 1];
    B_long = [params.Ts^2/2; params.Ts];
    A_lat = A_long;
    B_lat = B_long;
    params.A_discrete = blkdiag(A_long, A_lat);
    params.B_discrete = blkdiag(B_long, B_lat);
    
    % Initial condition
    params.egoInit = [0; 70/3.6; params.laneCenter2; 0];
    
    % Vehicle dimensions
    params.egoLength = 4;
    params.egoWidth = 2;
    
    %% Obstacles
    % Discrete dynamics
    A_long = [1 params.Ts; 0 1];
    B_long = [params.Ts^2/2; params.Ts];
    A_lat = A_long;
    B_lat = B_long;
    C_long = [1 0];
    C_lat = C_long;
    
    % Longitudinal velocity tracking
    costQ = 1;
    costR = 1;
    A_temp = A_long(2,2);
    B_temp = B_long(2);
    K_long = dlqr(A_temp, B_temp, costQ, costR);
    
    % Lateral position tracking
    costQ = diag([1 1]);
    costR = 1;
    K_lat = dlqr(A_lat, B_lat, costQ, costR);
    
    % Closed-loop dynamics
%     A_long_CL = [A_long(1,:); A_long(2,1) A_temp-B_temp*K_long];
    A_long_CL = A_long;
    A_lat_CL = A_lat - B_lat * K_lat;
    A_CL = blkdiag(A_long_CL, A_lat_CL);
%     B_long_CL = [zeros(1,2); 0 eye(1) - (A_temp - B_temp * K_long)];
    B_long_CL = zeros(2,2);
    B_lat_CL = eye(2) - A_lat_CL;
    B_CL = blkdiag(B_long_CL, B_lat_CL);
    C_CL = blkdiag(C_long, C_lat);
    
    % State-space dimensions
    obstacle.nx = 4;
    obstacle.nu = 2;
    obstacle.ny = 2;
    
    % State-space matrices
    obstacle.E = A_CL;
    obstacle.F = B_CL * [0; 70/3.6; params.laneCenter2; 0];
    obstacle.H = C_CL;
    
    % Noise
%     obstacle.G = B_CL * diag([0 (10/3.6/5)^2 (params.laneWidth/4/5)^2 0]) * B_CL';
    obstacle.G = blkdiag(diag([0 (1/5)^2]), B_lat_CL * diag([(params.laneWidth/4/5)^2 0]) * B_lat_CL');
    obstacle.initMean = [0; 70/3.6; params.laneCenter1; 0];
    obstacle.initVar = diag([0 (20/3.6/5)^2 0 (5/3.6/5)^2]);
    obstacle.measNoiseMean = [0; 0];
    obstacle.measNoiseVar = diag([(5/5)^2 (1/5)^2]);
    
    % Vehicle dimensions
    obstacle.Length = 4;
    obstacle.Width = 2;
    
    % Face coefficients
    obstacle.numFaces = 4;
    obstacle.aCoeff = [1 0; -1 0; 0 -1; 0 1];
    obstacle.cCoeff = [-1 0 0 0; 1 0 0 0; 0 0 1 0; 0 0 -1 0];
    obstacle.dCoeff = [-obstacle.Length/2; -obstacle.Length/2; -obstacle.Width/2; -obstacle.Width/2];
    
    % Data structure
    params.obstacles{1} = obstacle;

    %% Cost
    % Target
    params.targetState = [0; 70/3.6; params.laneCenter2; 0];
    
    % Weights
    params.targetStateWeight = diag([0 1 1 0]);
    params.inputWeight = diag([0.1 0.1]);

    %% Chance constraints
    % Safety margin
    params.epsilon = 0.05;
    
    % Certainty
    params.beta = 1e-3;
    
    % Number of samples
    params.Ns = 1e3;

    %% Other
    params.bigM = 1e4;
    
end