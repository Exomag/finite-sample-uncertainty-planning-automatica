function generatePlots()
    %% Simulation
    % RNG seed
    rng(0);
    
    % Load parameters
    params = loadParams();

    % Instantiate horizon data
    stepsHorizon = params.N;
    tHorizon = 0 : params.Ts : params.Ts * stepsHorizon;

    % Instantiate ego vehicle data
    xEgo = zeros(params.nx, stepsHorizon+1);
    uEgo = zeros(params.nu, stepsHorizon);

    % Instantiate ostacle data
    numObstacles = length(params.obstacles);
    xObs = cell(numObstacles, 1);
    xObsMeas = cell(numObstacles, 1);
    for iObs = 1 : length(xObs)
        currObs = params.obstacles{iObs};
        xObs{iObs} = zeros(currObs.nx, stepsHorizon+1);
        xObsMeas{iObs} = zeros(currObs.ny, stepsHorizon);
    end

    % Instantiate Kalman filter data
    xObsMeanPrior = cell(numObstacles, 1);
    xObsMeanPosterior = cell(numObstacles, 1);
    xObsVarPrior = cell(numObstacles, 1);
    xObsVarPosterior = cell(numObstacles, 1);
    for iObs = 1 : length(xObs)
        currObs = params.obstacles{iObs};
        xObsMeanPrior{iObs} = zeros(currObs.nx, stepsHorizon);
        xObsMeanPosterior{iObs} = zeros(currObs.nx, stepsHorizon+1);
        xObsVarPrior{iObs} = zeros(currObs.nx, currObs.nx, stepsHorizon);
        xObsVarPosterior{iObs} = zeros(currObs.nx, currObs.nx, stepsHorizon+1);
    end

    % Instantiate optimization data
    paramsData = cell(stepsHorizon, 1);
    problemData = cell(stepsHorizon, 1);
    solutionData = cell(stepsHorizon, 1);
    totalCost = 0;

    % Initialize ego vehicle state
    xEgo(:,1) = params.egoInit;

    % Initialize obstacle states
    for iObs = 1 : length(xObs)
        xObs{iObs}(:,1) = params.obstacles{iObs}.initMean;
    end

    % Initialize Kalman filter state
    for iObs = 1 : length(xObs)
        currObs = params.obstacles{iObs};
        xObsMeanPosterior{iObs}(:,1) = currObs.initMean;
        xObsVarPosterior{iObs}(:,:,1) = currObs.initVar;
    end

    % Receding Horizon
    for iHor = 1 : stepsHorizon
        % Update parameters
        params.egoInit = xEgo(:,iHor);
        for iObs = 1 : length(xObs)
            params.obstacles{iObs}.initMean = xObsMeanPosterior{iObs}(:,iHor);
            params.obstacles{iObs}.initVar = xObsVarPosterior{iObs}(:,:,iHor);
        end
        paramsData{iHor} = params;

        % Problem data
        problemData{iHor} = getProblemData(params);

        % Solve optimization problem
        solution = solveProblem(params, problemData{iHor});
        solutionData{iHor} = solution;

        % Update states
        uEgo(:,iHor) = solution.uVal(:,1);
        xEgo(:,iHor+1) = params.A_discrete * xEgo(:,iHor) + params.B_discrete * uEgo(:,iHor);
        for iObs = 1 : length(xObs)
            currObs = params.obstacles{iObs};
            w = mvnrnd(zeros(1, currObs.nx), currObs.G)';
            xObs{iObs}(:,iHor+1) = currObs.E * xObs{iObs}(:,iHor) + currObs.F + w;
        end

        % Total cost
        totalCost = totalCost + uEgo(:,iHor)' * params.inputWeight * uEgo(:,iHor) + ...
            (xEgo(:,iHor+1) - params.targetState)' * params.targetStateWeight * (xEgo(:,iHor+1) - params.targetState);

        % Measurements
        for iObs = 1 : length(xObs)
            currObs = params.obstacles{iObs};
            v = mvnrnd(currObs.measNoiseMean, currObs.measNoiseVar)';
            xObsMeas{iObs}(:,iHor) = currObs.H * xObs{iObs}(:,iHor+1) + v;
        end

        % Kalman filter
        for iObs = 1 : length(xObs)
            currObs = params.obstacles{iObs};
            xObsMeanPrior{iObs}(:,iHor) = currObs.E * xObsMeanPosterior{iObs}(:,iHor) + currObs.F;
            xObsVarPrior{iObs}(:,:,iHor) = currObs.E * xObsVarPosterior{iObs}(:,:,iHor) * currObs.E' + currObs.G;
            y = xObsMeas{iObs}(:,iHor) - currObs.H * xObsMeanPrior{iObs}(:,iHor);
            S = currObs.H * xObsVarPrior{iObs}(:,:,iHor) * currObs.H' + currObs.measNoiseVar;
            K = xObsVarPrior{iObs}(:,:,iHor) * currObs.H' * inv(S);
            xObsMeanPosterior{iObs}(:,iHor+1) = xObsMeanPrior{iObs}(:,iHor) + K * y;
            xObsVarPosterior{iObs}(:,:,iHor+1) = (eye(currObs.nx) - K * currObs.H) * xObsVarPrior{iObs}(:,:,iHor);
        end
    end

    %% Plots
    % Figure 1 - Frames
    figure();
    indFrames = 1 : 5 : params.N;
    numFrames = length(indFrames);
    obstacle = params.obstacles{1};
    xlimMin = -5;
    xlimMax = 85;
    for i = 1 : numFrames
        subplot(numFrames, 1, i); hold on; grid on;

        currInd = indFrames(i);
        plot(xEgo(1,1:currInd), xEgo(3,1:currInd), 'b-o');
        plot(solutionData{currInd}.xVal(1,:), solutionData{currInd}.xVal(3,:), '-s', 'Color', [0 100/255 0]);
        plot(xObs{1}(1,1:currInd), xObs{1}(3,1:currInd), 'r-');

        rectangle('Position', [xEgo(1,currInd)-params.egoLength/2, xEgo(3,currInd)-params.egoWidth/2, params.egoLength, params.egoWidth],'EdgeColor', 'b', 'FaceColor', 'none', 'LineWidth', 0.25)
        rectangle('Position', [xObs{1}(1,currInd)-obstacle.Length/2, xObs{1}(3,currInd)-obstacle.Width/2, obstacle.Length, obstacle.Width],'EdgeColor', 'r', 'FaceColor', 'none', 'LineWidth', 0.25)

        plot([xlimMin xlimMax], [params.laneWidth params.laneWidth], '-k', 'LineWidth', 2);
        plot([xlimMin xlimMax], [-params.laneWidth -params.laneWidth], 'k', 'LineWidth', 2);
        plot([xlimMin xlimMax], [0 0], '--k', 'LineWidth', 1);

        xlim([xlimMin xlimMax]);
        if i ~= numFrames
            set(gca,'xticklabel',[])
        else
            xlabel('$$x_1$$ (\si{\meter})', 'interpreter', 'latex');
        end

        ylim([-params.laneWidth-0.5 params.laneWidth+0.5]);
        ylabel('$$x_2$$ (\si{\meter})', 'interpreter', 'latex');

        if i == 1
            legend({'past', 'plan', 'adversary'}, 'Location', 'East', 'Orientation', 'horizontal', 'interpreter', 'latex');
        else
            legend('hide');
        end

        set(gca,'TickLabelInterpreter', 'latex')
    end

    save2tikz('plots/ClosedLoop_Frames')

    % Figure - Inputs
    figure();
    subplot(2, 1, 1); hold on; grid on;
    plot(tHorizon(1:end-1), uEgo(1,:), 'b-o');
    plot(tHorizon(1:end-1), solutionData{1}.uVal(1,:), '-s', 'Color', [0 100/255 0]);
    plot(tHorizon(1:end-1), params.xAccMin * ones(size(tHorizon(1:end-1))), '--r');
    plot(tHorizon(1:end-1), params.xAccMax * ones(size(tHorizon(1:end-1))), '--r');
    xlim([tHorizon(1) tHorizon(end-1)]);
    ylim([params.xAccMin params.xAccMax]); ylabel('$$u_1$$ (\si{\meter\per\square\second})', 'interpreter', 'latex');
    legend({'closed loop', 'open loop'}, 'Location', 'NorthWest', 'interpreter', 'latex');
    set(gca,'TickLabelInterpreter', 'latex');
    subplot(2, 1, 2); hold on; grid on;
    plot(tHorizon(1:end-1), uEgo(2,:), 'b-o');
    plot(tHorizon(1:end-1), solutionData{1}.uVal(2,:), '-s', 'Color', [0 100/255 0]);
    plot(tHorizon(1:end-1), params.yAccMin * ones(size(tHorizon(1:end-1))), '--r');
    plot(tHorizon(1:end-1), params.yAccMax * ones(size(tHorizon(1:end-1))), '--r');
    xlim([tHorizon(1) tHorizon(end-1)]); xlabel('Time $$t$$ (\si{\second})', 'interpreter', 'latex');
    ylim([params.yAccMin params.yAccMax]); ylabel('$$u_2$$ (\si{\meter\per\square\second})', 'interpreter', 'latex');
    legend('hide');
    set(gca,'TickLabelInterpreter', 'latex');

    save2tikz('plots/ClosedLoop_Inputs')

    %% Figure - Obstacle position
    figure();
    meanOpenLoop = mean(squeeze(problemData{1}.obstacleStateSamples{1}(2,:,:)));
    stdOpenLoop = sqrt(var(squeeze(problemData{1}.obstacleStateSamples{1}(2,:,:))));
    meanClosedLoop = xObsMeanPosterior{1}(2,:);
    stdClosedLoop = sqrt(squeeze(xObsVarPosterior{1}(2,2,:)))';
    subplot(2, 1, 1); hold on; grid on;
    plot(tHorizon, meanClosedLoop, 'b-o');
    plot(tHorizon, meanOpenLoop, '-s', 'Color', [0 100/255 0]);
    plot(tHorizon, meanClosedLoop + stdClosedLoop, 'b--');
    plot(tHorizon, meanClosedLoop - stdClosedLoop, 'b--');
    plot(tHorizon, meanOpenLoop + stdOpenLoop, '--', 'Color', [0 100/255 0]);
    plot(tHorizon, meanOpenLoop - stdOpenLoop, '--', 'Color', [0 100/255 0]);
    xlim([tHorizon(1) tHorizon(end-1)]);
    ylabel('$$\chi_3$$ (\si{\meter\per\second})', 'interpreter', 'latex');
    legend({'closed loop', 'open loop'}, 'Location', 'SouthEast', 'interpreter', 'latex');
    set(gca,'TickLabelInterpreter', 'latex');
    meanOpenLoop = mean(squeeze(problemData{1}.obstacleStateSamples{1}(4,:,:)));
    stdOpenLoop = sqrt(var(squeeze(problemData{1}.obstacleStateSamples{1}(4,:,:))));
    meanClosedLoop = xObsMeanPosterior{1}(4,:);
    stdClosedLoop = sqrt(squeeze(xObsVarPosterior{1}(4,4,:)))';
    subplot(2, 1, 2); hold on; grid on;
    plot(tHorizon, meanClosedLoop, 'b-o');
    plot(tHorizon, meanOpenLoop, '-s', 'Color', [0 100/255 0]);
    plot(tHorizon, meanClosedLoop + stdClosedLoop, 'b--');
    plot(tHorizon, meanClosedLoop - stdClosedLoop, 'b--');
    plot(tHorizon, meanOpenLoop + stdOpenLoop, '--', 'Color', [0 100/255 0]);
    plot(tHorizon, meanOpenLoop - stdOpenLoop, '--', 'Color', [0 100/255 0]);
    xlim([tHorizon(1) tHorizon(end-1)]); xlabel('Time $$t$$ (\si{\second})', 'interpreter', 'latex');
    ylabel('$$\chi_4$$ (\si{\meter\per\second})', 'interpreter', 'latex');
    legend('hide');
    set(gca,'TickLabelInterpreter', 'latex');
    save2tikz('plots/ClosedLoop_Uncertainty')
end
